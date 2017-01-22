% DEMJE2 Demonstrate nonlinear time series joint estimation for Mackey-Glass chaotic time series
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
%   GSSM_MACKEY_GLASS
%   Copyright (c) Oregon Health & Science University (2006)
%
%   This file is part of the ReBEL Toolkit. The ReBEL Toolkit is available free for
%   academic use only (see included license file) and can be obtained from
%   http://choosh.csee.ogi.edu/rebel/.  Businesses wishing to obtain a copy of the
%   software should contact rebel@csee.ogi.edu for commercial licensing information.
%
%   See LICENSE (which should be part of the main toolkit distribution) for more
%   detail.

%===============================================================================================

clc;
clear all; close all;

fprintf('\nDEMJE2:  This demonstration shows how the ReBEL toolkit is used for joint estimation\n');
fprintf('         on a nonlinear time series (Mackey-Glass-30) problem. The scalar observation\n');
fprintf('         is corrupted by additive white Gaussian noise. A neural network is used as a\n');
fprintf('         generative model for the time series. We estimate both the model parameters and\n');
fprintf('         the underlying clean state from the noisy observations.\n');
fprintf('         We compare the performance of an EKF and a SRCDKF by iterating on the same sequence.\n\n');
fprintf('    NOTE : This demos is quite computationally expensive... so on a slow computer it might take a while.\n\n');


%--- General setup

addrelpath('../gssm');         % add relative search path to example GSSM files to MATLABPATH
addrelpath('../data');         % add relative search path to example data files to MATLABPATH

%--- Initialise GSSM model from external system description script.
model = gssm_mackey_glass('init');


%--- Load normalized Mackey glass data set

load('mg30_normalized.mat');                           % loads mg30_data from ../data/mg30_normalized.mat

mg30_data = mg30_data(100:100+300-1);


%--- Build state space data matrix of input data

X = datamat(mg30_data, model.statedim);                 % pack vector of data into datamtrix for NN input

[dim,N]  = size(X);                                     % dimension and number of datapoints
y  = zeros(model.obsdim,N);                             % observation data buffer

clean_signal_var = var(mg30_data);                      % determine variance of clean time series

SNR = 3;                                                % 3db SNR
onoise_var = clean_signal_var/10^(SNR/10);              % determine needed observation noise variance for a given SNR

model.oNoise.cov = onoise_var;                            % set observation noise covariance

onoise = model.oNoise.sample( model.oNoise, N);   % generate observation noise

y   = model.hfun( model, X, onoise);   % generate observed time series (corrupted with observation noise)

%----

ftype1 = 'ekf';
ftype2 = 'srcdkf';


%--- Setup argument data structure which serves as input to
%--- the 'geninfds' function. This function generates the InferenceDS and
%--- SystemNoiseDS data structures which are needed by all inference algorithms
%--- in the PiLab toolkit.

Arg.type = 'joint';                                  % inference type (state estimation)
Arg.tag = 'Joint estimation for GSSM_MACKEY_GLASS system.';  % arbitrary ID tag
Arg.model = model;                                   % GSSM data structure of external system

InfDS = geninfds(Arg);                               % create inference data structure

[pNoise1, oNoise1, InfDS1] = gensysnoiseds(InfDS,ftype1);    % generate process and observation noise sources for EKF
[pNoise2, oNoise2, InfDS2] = gensysnoiseds(InfDS,ftype2);    % generate process and observation noise sources for SRCDKF


%--- Setup runtime buffers

Xh = zeros(InfDS.statedim,N);          % state estimation buffer
Px = eye(InfDS.statedim);            % initial state covariance
Px(model.statedim+1:end,model.statedim+1:end) = 0.1*eye(model.paramdim);

Xh(model.statedim+1:end,1) = mlpweightinit(model.nodes);              % randomize initial model parameters

Xh1 = Xh;
Px1 = Px;

Xh2 = Xh;
Sx2 = chol(Px)';                     % SRCDKF is a square-root algorithm and hence it operates on the Cholesky factor
                                     % of the covariance matrix
number_of_runs = 10;                  % we will iterate over the data 'number_of_runs' times

mse1 = zeros(1,number_of_runs);       % buffers to store the MSE of each runs estimate
mse2 = mse1;

mse1(1) = mean((y(1,:)-X(1,:)).^2)/var(y(1,:));    % initial MSE of noisy signal
mse2(1) = mse1(1);

%--- Setup process noise data structures for joint estimation

  pNoiseAdaptMethod = 'anneal';                                % setup process noise adaptation method (improves convergence)
  pNoiseAdaptParams = [0.995 1e-7];                            % annealing factor = 0.95     annealing floor variance = 1e-8

  pNoiseCov0 = 1e-4*eye(model.paramdim);

  pNoise1.adaptMethod = pNoiseAdaptMethod;
  pNoise1.adaptParams = pNoiseAdaptParams;

  pNoise2.adaptMethod = pNoiseAdaptMethod;
  pNoise2.adaptParams = pNoiseAdaptParams;

  pNoise1.cov(2:end,2:end) = pNoiseCov0;         % set initial variance of process noise parameter estimation subvector
  pNoise2.cov(2:end,2:end) = chol(pNoiseCov0)';  % set initial variance of process noise parameter estimation subvector


%---

fprintf('\n Running joint estimators ... \n\n');


%--- Call inference algorithm / estimator

for k=1:number_of_runs,

  fprintf(' [%d:%d] ',k,number_of_runs);


  %------------------- Extended Kalman Filter ------------------------------------


  [Xh1, Px1, pNoise1] = ekf(Xh1(:,1), Px1, pNoise1, oNoise1, y, [], [], InfDS1);


  %------------------- Square-root Central Difference Kalman Filter -------------

  InfDS2.spkfParams = sqrt(3); ;                                 % scale factor (CDKF parameter)

  [Xh2, Sx2, pNoise2] = srcdkf(Xh2(:,1), Sx2, pNoise2, oNoise2, y, [], [], InfDS2);

  %---------------------------------------------------------------------------------


  %--- Calculate normalized mean square estimation error

  mse1(k+1) = mean((Xh1(1,:)-X(1,:)).^2)/var(y(1,:));
  mse2(k+1) = mean((Xh2(1,:)-X(1,:)).^2)/var(y(1,:));

  %--- Plot results

  figure(1); clf; subplot('position',[0.025 0.1 0.95 0.8]);
  p1 = plot(X(1,:),'b','linewidth',2); hold on
  p2 = plot(y,'g+');
  p3 = plot(Xh1(1,:),'m');
  p4 = plot(Xh2(1,:),'r'); hold off
  legend([p1 p2 p3 p4],'clean','noisy','EKF estimate','SRCDKF estimate');
  xlabel('time');
  ylabel('x');
  title('DEMSE3 : Mackey-Glass-30 Chaotic Time Series Joint Estimation');

  figure(2);
  p1 = plot(mse1(2:k+1),'m-o'); hold on;
  p2 = plot(mse2(2:k+1),'r-s'); hold off;
  legend([p1 p2],'EKF','SRCDKF');
  title('Normalized MSE of Estimates');
  xlabel('k');
  ylabel('MSE');
  drawnow

  fprintf('  Mean-square-error (MSE) of estimates : EKF = %4.3f    SRCDKF = %4.3f\n', mse1(k+1), mse2(k+1));


  %-- Copy last estimate of model parameters to initial buffer position for next iteration...

  Xh1(model.statedim+1:end,1) = Xh1(model.statedim+1:end,end);              % copy model parameters over
  Xh1(1:model.statedim,1) = zeros(model.statedim,1);                        % reset state estimate
  Px1_temp = eye(InfDS.statedim);                                           % copy covariance of parameter estimates
  Px1_temp(model.statedim+1:end,model.statedim+1:end) = Px1(model.statedim+1:end,model.statedim+1:end);
  Px1 = Px1_temp;

  Xh2(model.statedim+1:end,1) = Xh2(model.statedim+1:end,end);              % copy model parameters over
  Xh2(1:model.statedim,1) = zeros(model.statedim,1);                        % reset state estimate
  Sx2_temp = eye(InfDS.statedim);                                           % copy covariance of parameter estimates
  Sx2_temp(model.statedim+1:end,model.statedim+1:end) = Sx2(model.statedim+1:end,model.statedim+1:end);
  Sx2 = Sx2_temp;


end


%--- House keeping

remrelpath('../gssm');       % remove relative search path to example GSSM files from MATLABPATH
remrelpath('../data');       % remove relative search path to example data files from MATLABPATH
