% DEMSE3  Demonstrate nonlinear time series state estimation for Mackey-Glass chaotic time series
%
%  The Mackey-Glass time-delay differential equation is defined by
%
%            dx(t)/dt = 0.2x(t-tau)/(1+x(t-tau)^10) - 0.1x(t)
%
%  When x(0) = 1.2 and tau = 17, we have a non-periodic and non-convergent time series that
%  is very sensitive to initial conditions. (We assume x(t) = 0 when t < 0.)
%
%  We assume that the chaotic time series is generated with by a nonlinear autoregressive
%  model where the nonlinear functional unit is a feedforward neural network. We use a
%  tap length of 6 and a 6-4-1 MLP neural network (using the Netlab toolkit) with hyperbolic
%  tangent activation functions in the hidden layer and a linear output activation.
%
%   See also
%   GSSM_MACKEY_GLASS, DEMSE1, DEMSE2
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

fprintf('\nDEMSE3 : Demonstrate nonlinear state estimation for Mackey-Glass chaotic time series\n\n');


%--- General setup

addrelpath('../gssm');         % add relative search path to example GSSM files to MATLABPATH
addrelpath('../data');         % add relative search path to example data files to MATLABPATH

%--- Initialise GSSM model from external system description script.

model = gssm_mackey_glass('init');

%--- Load normalized Mackey glass data set

load('mg30_normalized.mat');                            % load 'mg30_data' variable

mg30_data = mg30_data(1:1000);                          % only use 1000 data points


%--- Build state space data matrix of input data

X = datamat(mg30_data, model.statedim);                 % pack vector of data into datamtrix for NN input

[dim,N]  = size(X);                                     % dimension and number of datapoints
y  = zeros(model.obsdim,N);                             % observation data buffer

clean_signal_var = var(mg30_data);                      % determine variance of clean time series

SNR = 3;                                                % 3db SNR
onoise_var = clean_signal_var/10^(SNR/10);              % determine needed observation noise variance for a given SNR

model.oNoise.cov = onoise_var;                            % set observation noise covariance

onoise = model.oNoise.sample( model.oNoise, N);   % generate observation noise

y   = model.hfun( model, X, onoise);    % generate observed time series (corrupted with observation noise)

figure(1);
p1=plot(X(1,:),'b'); hold on;
p2=plot(y,'g+');
legend([p1 p2],'clean','noisy');
xlabel('time - k');
drawnow

%--- Ask the user which inference algorithm to use
% ftype = input('Type of estimator [ ekf, ukf, cdkf, srcdkf or srukf ] ? ','s');

%--- Use a couple of different filters...
lftype = {'ekf','ukf','cdkf','srukf','srcdkf'};

for k=1:5,
    
  ftype = lftype{k};    

  %--- Setup argument data structure which serves as input to
  %--- the 'geninfds' function. This function generates the InferenceDS and
  %--- SystemNoiseDS data structures which are needed by all inference algorithms
  %--- in the PiLab toolkit.

  Arg.type = 'state';                                  % inference type (state estimation)
  Arg.tag = 'State estimation for GSSM_MACKEY_GLASS system.';  % arbitrary ID tag
  Arg.model = model;                                   % GSSM data structure of external system

  InfDS = geninfds(Arg);                               % Create inference data structure and

  [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS,ftype);       % generate process and observation noise sources


  %--- Setup runtime buffers

  Xh = zeros(InfDS.statedim,N);          % state estimation buffer
  Xh(:,1) = X(:,1);     % initial estimate of state E[X(0)]
  Px = eye(InfDS.statedim);              % initial state covariance


  %--- Call inference algorithm / estimator

  switch ftype


    %------------------- Extended Kalman Filter ------------------------------------
    case 'ekf'

        [Xh, Px] = ekf(Xh(:,1), Px, pNoise, oNoise, y, [], [], InfDS);


    %------------------- Unscented Kalman Filter -----------------------------------
    case 'ukf'

        alpha = 1;         % scale factor (UKF parameter)
        beta  = 2;         % optimal setting for Gaussian priors (UKF parameter)
        kappa = 0;         % optimal for state dimension=2 (UKF parameter)

        InfDS.spkfParams = [alpha beta kappa];

        [Xh, Px] = ukf(Xh(:,1), Px, pNoise, oNoise, y, [], [], InfDS);


    %------------------- Central Difference Kalman Filter ---------------------------
    case 'cdkf'

        InfDS.spkfParams = sqrt(3);    % scale factor (CDKF parameter h)

        [Xh, Px] = cdkf(Xh(:,1), Px, pNoise, oNoise, y, [], [], InfDS);


    %------------------- Square Root Unscented Kalman Filter ------------------------
    case 'srukf'

        alpha = 1;         % scale factor (UKF parameter)
        beta  = 2;         % optimal setting for Gaussian priors (UKF parameter)
        kappa = 0;         % optimal for state dimension=2 (UKF parameter)

        Sx = chol(Px)';

        InfDS.spkfParams = [alpha beta kappa];

        [Xh, Sx] = srukf(Xh(:,1), Sx, pNoise, oNoise, y, [], [], InfDS);


    %------------------- Square Root Central Difference Kalman Filter ---------------
    case 'srcdkf'

        InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)

        Sx = chol(Px)';

        [Xh, Sx] = srcdkf(Xh(:,1), Sx, pNoise, oNoise, y, [], [], InfDS);


   otherwise

    error(' Unknown estimator!');

  end

  %--- Plot results

  figure(k); clf;
  p1 = plot(X(1,:)); hold on
  p2 = plot(y,'g+');
  p3 = plot(Xh(1,:),'r'); hold off;
  legend([p1 p2 p3],'clean','noisy',[ftype ' estimate']);
  xlabel('time');
  title('DEMSE3 : Mackey-Glass-30 Chaotic Time Series State Estimation');


  %--- Calculate mean square estimation error

  mse = mean((Xh(1,:)-X(1,:)).^2);
  disp([ftype ' : Mean square error (MSE) of estimate : ' num2str(mse)]);

end
  
%--- House keeping

remrelpath('../gssm');       % remove relative search path to example GSSM files from MATLABPATH
remrelpath('../data');       % remove relative search path to example data files from MATLABPATH