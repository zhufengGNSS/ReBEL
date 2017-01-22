% DEMPE1  Demonstrate parameter estimation on a simple 2nd order LTI system.
%
%   This is a simple demonstration of how to use the ReBEL toolkit for parameter estimation on
%   a simple 2nd order LTI system.
%
%   See also
%   DEMSE1, DEMJE1, GSSM_LTI1
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

clc;
clear all;

fprintf('\nDEMPE1 : This demonstration shows how the ReBEL toolkit is used for parameter estimation\n');
fprintf('         on a 2nd order LTI system. The scalar observation is corrupted by additive white\n');
fprintf('         Gaussian noise.\n\n');


%--- General setup

addrelpath('../gssm');         % add relative search path to example GSSM files to MATLABPATH
addrelpath('../data');         % add relative search path to example data files to MATLABPATH

%--- Initialise GSSM model from external system description script.

model = gssm_lti1('init');


%--- Generate some data

N  = 600;                                               % number of datapoints
X  = zeros(model.statedim,N);                           % state data buffer
y  = zeros(model.obsdim,N);                             % observation data buffer

pnoise = model.pNoise.sample( model.pNoise, N);   % generate process noise
onoise = model.oNoise.sample( model.oNoise, N);   % generate observation noise

X(:,1) = [1 0]';                                          % initial state
y(1)   = model.hfun( model, X(:,1), onoise(1), []); % observation of initial state
for j=2:N,
  X(:,j) = model.ffun( model, X(:,j-1), pnoise(j-1), []);
  y(j)   = model.hfun( model, X(:,j), onoise(j), []);
end


%--- Ask the user which inference algorithm to use

ftype = input('Type of estimator [ kf / ekf / ukf / cdkf / srukf / srcdkf ] ? ','s');

if ~stringmatch(ftype,{'kf','ekf','ukf','cdkf','srukf','srcdkf'})
    error('Unknown estimator!');
end


%--- Setup argument data structure which serves as input to
%--- the 'geninfds' function. This function generates the InferenceDS and
%--- SystemNoiseDS data structures which are needed by all inference algorithms
%--- in the PiLab toolkit.

Arg.type = 'parameter';                              % inference type (parameter estimation)
Arg.tag = 'Parameter estimation for GSSM_LTI1 system.';  % arbitrary ID tag
Arg.model = model;                                   % GSSM data structure of external system
Arg.paramFunSelect = 'ffun';                         % Only use gssm_lti's FFUN for parameter estimation purposes
Arg.paramFFunOutIdxVec = [1];                        % Only use the first element of the output vector of FFUN as the
                                                     %  observation. This is for auto-regressive state space models.

InfDS = geninfds(Arg);                                % Create inference data structure
[pNoise, oNoise, InfDS] = gensysnoiseds(InfDS,ftype); % Generate process and observation noise data structures

%--- Setup runtime buffers

Wh = zeros(InfDS.statedim,N-1);                      % meta system state estimation buffer
Wh(:,1) = 0.01*randn(InfDS.statedim,1);              % initial condition of state E[Wh(0)]
Pw = eye(InfDS.statedim);                            % initial state covariance
U2 = X(:,1:end-1);                                   % exogenous input for observation function of meta system
Y  = X(InfDS.paramFFunOutIdxVec,2:end);              % meta system observations


%--- Call inference algorithm / estimator

switch ftype

    %------------------- Linear Kalman Filter --------------------------------------
    case 'kf'

        pNoise.adaptMethod = 'anneal';             % process noise adaptation method : annealing
        pNoise.adaptParams = [0.9 1e-8];           % annealing parameters [anneal_factor minimum_variance]

        [Wh, Pw, pNoise] = kf(Wh(:,1), Pw, pNoise, oNoise, Y, [], U2, InfDS);

    %------------------- Extended Kalman Filter ------------------------------------
    case 'ekf'

        pNoise.adaptMethod = 'anneal';
        pNoise.adaptParams = [0.9 1e-8];

        [Wh, Pw, pNoise] = ekf(Wh(:,1), Pw, pNoise, oNoise, Y, [], U2, InfDS);


    %------------------- Square Root Unscented Kalman Filter -----------------------------------
    case 'ukf'

        alpha = 1;        % scale factor (UKF parameter)
        beta  = 3;        % optimal setting for Gaussian priors (UKF parameter)
        kappa = 0;        % optimal for state dimension=2 (UKF parameter)

        InfDS.spkfParams = [alpha beta kappa]; % UKF parameter vector

        pNoise.adaptMethod = 'anneal';
        pNoise.adaptParams = [0.9 1e-8];

        [Wh, Pw, pNoise] = ukf(Wh(:,1), Pw, pNoise, oNoise, Y, [], U2, InfDS);


    %------------------- Square Root Central Difference Kalman Filter ---------------------------
    case 'cdkf'

        InfDS.spkfParams = sqrt(3);     % scale factor (CDKF parameter h)

        pNoise.adaptMethod = 'anneal';
        pNoise.adaptParams = [0.9 1e-8];


        [Wh, Pw, pNoise] = cdkf(Wh(:,1), Pw, pNoise, oNoise, Y, [], U2, InfDS);


    %------------------- Square Root Unscented Kalman Filter -----------------------------------
    case 'srukf'

        alpha = 1;        % scale factor (UKF parameter)
        beta  = 3;        % optimal setting for Gaussian priors (UKF parameter)
        kappa = 0;        % optimal for state dimension=2 (UKF parameter)

        InfDS.spkfParams = [alpha beta kappa]; % UKF parameter vector

        pNoise.adaptMethod = 'anneal';
        pNoise.adaptParams = [0.9 1e-8];

        Sw = chol(Pw)';

        [Wh, Sw, pNoise] = srukf(Wh(:,1), Sw, pNoise, oNoise, Y, [], U2, InfDS);


    %------------------- Square Root Central Difference Kalman Filter ---------------------------
    case 'srcdkf'

        InfDS.spkfParams = sqrt(3);     % scale factor (CDKF parameter h)

        pNoise.adaptMethod = 'anneal';
        pNoise.adaptParams = [0.9 1e-8];

        Sw = chol(Pw)';

        [Wh, Sw, pNoise] = srcdkf(Wh(:,1), Sw, pNoise, oNoise, Y, [], U2, InfDS);

    %---
    otherwise error('That estimator/filter type is not recognized.');

end


%--- Display results

true_model_trace = cvecrep(model.params,N-1);

figure(1); clf;
p11=plot(true_model_trace(1,:),'b--','linewidth',3); hold on;
p12=plot(true_model_trace(2,:),'r--','linewidth',3);
p21=plot(Wh(1,:),'b');
p22=plot(Wh(2,:),'r'); hold off;
axis([1 N -1.5 2.5]);
legend([p11 p21 p12 p22],'true parameter 1','estimated parameter 1','true parameter 2','estimated parameter 2',0);
ylabel('model parameters')
xlabel('time');
title(['DEMPE1 : LTI System Parameter Estimation  ( ' ftype ' )']);


%--- House keeping

remrelpath('../gssm');       % remove relative search path to example GSSM files from MATLABPATH
remrelpath('../data');       % remove relative search path to example data files from MATLABPATH
