%  DEMPE2  Demonstrate how the ReBEL toolkit is used to train a neural network
%          in efficient "low-level" mode by directly defining an InferenceDS
%          data structure.
%
%   We train a small 2 layer MLP on the standard XOR classification problem using
%   a Square-Root Central Difference Kalman Filter (SPKF variant)
%
%   --------------------------- NOTE ----------------------------------------------------
%
%   In this example we make use of a directly defined InferenceDS data structure, without
%   the use of a GSSM (general state space model) file. By doing this we are removing one
%   level of abstraction from the inference problem, thereby gaining execution speed. ReBEL
%   in general makes use of a two-layer abstraction approach to seperate the problem
%   definition (gssm file) from the actual inference/estimation algorithms. This allows for
%   easy state-,parameter and/or joint estimation to be done on the same model without having
%   to modify the underlying state-space formulation and/or estimator implementations. This
%   increased generality and ease of implementation however comes at the cost of increased
%   computational overhead, since the estimation algorithms access the underlying model functions
%   (as defined in gssm) via the InferenceDS abstraction/state-space-mapping layers.
%   This combined with Matlab's inherent function-calling overhead causes a less than desireable
%   speed penalty to be paid. But, c'est la vie... what we loose in speed we gain in protoyping
%   ease.
%
%   HOWEVER, once is free to describe your system directly in the InferenceDS layer (which is
%   what we'll be doing in this example) in order to lessen the function call overhead. The down
%   side is that you have to now make sure you implement the correct state-space reformulation
%   depending if you are doing state-, parameter or joint estimation and you also need to comply
%   with the interface expected by the different estimation algorithms.
%
%   This is thus an example showing how ReBEL can be used on a "lower level".
%
%   ---------------------------------------------------------------------------------------
%
%   See also
%
%   INFDS_TRAIN_NN
%
%   Copyright (c) Oregon Health & Science University (2006)
%
%   This file is part of the ReBEL Toolkit. The ReBEL Toolkit is available free for
%   academic use only (see included license file) and can be obtained from
%   http://choosh.csee.ogi.edu/rebel/.  Businesses wishing to obtain a copy of the
%   software should contact rebel@csee.ogi.edu for commercial licensing information.
%
%   See LICENSE (which should be part of the main toolkit distribution) for more
%   detail.

%=============================================================================================

clear all; clc;

fprintf('\nDEMPE2 : Demonstrate how the ReBEL toolkit is used to train a neural network\n');
fprintf('         in efficient "low-level" mode by directly defining an InferenceDS\n');
fprintf('         data structure. We train a small 2 layer MLP on the standard XOR\n');
fprintf('         classification problem using a Square-Root Central Difference Kalman\n');
fprintf('         Filter (SPKF variant). We do a single pass through the data set. This\n');
fprintf('         problem is not optimsed for optimal training and generalization perfor-\n');
fprintf('         mance. It is simply used to demonstrate the method.\n\n');


%--- General setup

addrelpath('../gssm');         % add relative search path to example GSSM files to MATLABPATH
addrelpath('../data');         % add relative search path to example data files to MATLABPATH

%--- Generate some data for the XOR problem

N=4000;                            % number of data points
Ntrain = 1000;                     % size of training set
Ntest  = N-Ntrain;                 % size of testing set

X = 1-2*rand(2,N);                 % uniformly draw input points from [-1 1 -1 1]
C = 0.95*sign(X(1,:).*X(2,:));      % standard XOR classes

C1Idx = find(C>0);
C2Idx = find(C<0);

figure(1);
clf;
plot(X(1,C1Idx),X(2,C1Idx),'b.'); hold on;
plot(X(1,C2Idx),X(2,C2Idx),'r.'); hold off;
xlabel('x1'); ylabel('x2');
title('XOR Input Data');
drawnow

Xtrain = X(:,1:Ntrain);              % trainet set
Ctrain = C(:,1:Ntrain);

Xtest = X(:,Ntrain+1:end);           % testing set
Ctest = C(:,Ntrain+1:end);


%--- Setup neural network and InferenceDS structures

InfDS = infds_train_nn;         % Directly generate InferenceDS data structure. See 'infds_train_nn.m' for detail

InfDS = fixinfds(InfDS);        % Make sure all required fields of InfDS are set and add required default fields which
                                % the user did not specify.

%--- Generate process and observation noise sources needed by the SPKF inference algorithm (srcdkf)
%    Since the SRCDKF is a square-root algorithm, all noise sources should be of type 'gaussian' and
%    cov_type 'sqrt'.

Arg.type = 'gaussian';                     % Gaussian process noise noise source
Arg.cov_type = 'sqrt';                     % Square root form
Arg.dim  = InfDS.statedim;                 % noise vector dimension
Arg.mu   = zeros(InfDS.statedim,1);        % zero mean
Arg.cov  = sqrt(1e-1)*eye(InfDS.statedim); % initial noise covariance (Cholesky factor used for square-root forms). Usually
                                           % a good idea to not set this to large initially. We will also aneal this during
                                           % training.

pNoise = gennoiseds(Arg);                 % generate process noise source with call to 'gennoiseds'

pNoise.adaptMethod = 'anneal';            % we will use the 'annealing' method to adapt the process noise covariance
pNoise.adaptParams = [0.95 1e-7];         % annealing factor = 0.98  and minimum allowed variance (variance floor) = 1e-8


Arg.type = 'gaussian';                    % Gaussian observation noise noise source
Arg.cov_type = 'full';                    % Set covariance matrix type : full covariance
Arg.dim  = InfDS.obsdim;                  % noise vector dimension
Arg.mu   = zeros(InfDS.obsdim,1);         % zero mean
Arg.cov  = eye(InfDS.obsdim);             % initial noise covariance (Cholesky factor used for square-root forms).
                                          % For parameter estimation the absolute value of the observation noise covariance
                                          % is not crucial. Only the relative values (across the outputs) determine the
                                          % relative weighting of the output errors.

oNoise = gennoiseds(Arg);                 % generate observation noise source with call to 'gennoiseds'


%--- Setup estimation buffers

Wh = zeros(InfDS.statedim, Ntrain);  % setup state buffer  (the NN parameters are the states of our state-space system
                                     % defined in 'infds_train_bb.m'


Wh(:,1) = mlpweightinit(InfDS.nodes);   % initialize initial parameter vector

Sw = eye(InfDS.statedim);               % Initial state covariance Cholesky factor


%--- Call estimator
%    Here we are calling the estimator in batch mode, i.e. we are passing it all the data (all observations) at once. The
%    estimator will recursively run through all the observations internally and return a vector (or matrix) of all the
%    estimates from k=1:N. This is a more efficient way of calling the estimator if we have all of the data available ofline.
%    The estimator can however also be called using a single observation per time instance (i.e. external recursion). This
%    would be the standard way of using the estimator in an on-line situation.

InfDS.spkfParams = sqrt(3);            % SPKF parameter : CDKF step size


%--- Call the SRCDKF estimator

[Wh, Sw, pNoise] = srcdkf(Wh(:,1), Sw, pNoise, oNoise, Ctrain, [], Xtrain, InfDS);  % train on the training set



%--- Calculate performance of trained neural network

NNparams = Wh(:,end);

Ytrain = mlpff(InfDS.olType, InfDS.nodes, Xtrain, NNparams);       % output on training set
Ytest  = mlpff(InfDS.olType, InfDS.nodes, Xtest, NNparams);        % output on testing set

%--- Plot classification results

figure(2); clf

subplot(211);
Y1Idx = find(Ytrain>0);
Y2Idx = find(Ytrain<0);
plot(Xtrain(1,Y1Idx),Xtrain(2,Y1Idx),'b.'); hold on;
plot(Xtrain(1,Y2Idx),Xtrain(2,Y2Idx),'r.'); hold off;
xlabel('x1'); ylabel('x2');
title('Classification Results on Training Set');

subplot(212);
Y1Idx = find(Ytest>0);
Y2Idx = find(Ytest<0);
plot(Xtest(1,Y1Idx),Xtest(2,Y1Idx),'b.'); hold on;
plot(Xtest(1,Y2Idx),Xtest(2,Y2Idx),'r.'); hold off;
xlabel('x1'); ylabel('x2');
title('Classification Results on Testing Set');

drawnow

%--- Calculate classification performance

cerror_train = sum(0.5*abs(sign(Ctrain)-sign(Ytrain)))/length(Ctrain);
cerror_test = sum(0.5*abs(sign(Ctest)-sign(Ytest)))/length(Ctest);

disp(['Classification error on training set : ' num2str(round(cerror_train*100)) ' %']);
disp(['Classification error on test set     : ' num2str(round(cerror_test*100)) ' %']);
disp(' ');


%--- House keeping

remrelpath('../gssm');       % remove relative search path to example GSSM files from MATLABPATH
remrelpath('../data');       % remove relative search path to example data files from MATLABPATH
