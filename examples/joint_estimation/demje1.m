% DEMJE1  Demonstrate joint estimation on a 2nd order LTI system.
%
%   This is a demonstration of how to use the ReBEL toolkit for joint estimation on
%   a simple 2nd order LTI system.
%
%   See also
%   GSSM_LTI1, DEMSE1, DEMPE1
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

%===============================================================================================

clc;
clear all; close all;

fprintf('\nDEMJE1 : This demonstration shows how the ReBEL toolkit is used for joint estimation\n');
fprintf('         on a 2nd order LTI system. The scalar observation is corrupted by additive white\n');
fprintf('         Gaussian noise. We compare the performance of the EKF to that of the SRUKF on the\n');
fprintf('         same sequence.\n\n');
fprintf('         Note : This example has not been fine tuned for optimal performance and can sometimes\n');
fprintf('                fail to convergence to a good solution. In such a case, simply re-run the experiment \n\n');


%--- General setup

addrelpath('../gssm');         % add relative search path to example GSSM files to MATLABPATH
addrelpath('../data');         % add relative search path to example data files to MATLABPATH

%--- Initialise GSSM model from external system description script.
model = gssm_lti1('init');


%--- Generate some data

N  = 1000;                                               % number of datapoints
X  = zeros(model.statedim,N);                           % state data buffer
y  = zeros(model.obsdim,N);                             % observation data buffer

pnoise = model.pNoise.sample( model.pNoise, N);   % generate process noise
onoise = model.oNoise.sample( model.oNoise, N);   % generate observation noise


X(:,1) = [1 0]';                                          % initial state
y(1)   = model.hfun( model, X(:,1), onoise(1), []); % observation of initial state
for j=2:N,
    X(:,j) = model.ffun( model, X(:,j-1), pnoise(:,j-1), []);
    y(j)   = model.hfun( model, X(:,j), onoise(:,j), []);
end


ftype1 = 'ekf';
ftype2 = 'srukf';


%--- Setup argument data structure which serves as input to
%--- the 'geninfds' function. This function generates the InferenceDS and
%--- SystemNoiseDS data structures which are needed by all inference algorithms
%--- in the PiLab toolkit.

Arg.type = 'joint';                                  % inference type (state estimation)
Arg.tag = 'Joint estimation for GSSM_LTI1 system.';  % arbitrary ID tag
Arg.model = model;                                   % GSSM data structure of external system

InfDS = geninfds(Arg);                               % create inference data structure

[pNoise1, oNoise1, InfDS1] = gensysnoiseds(InfDS,ftype1);    % generate process and observation noise sources for SRUKF
[pNoise2, oNoise2, InfDS2] = gensysnoiseds(InfDS,ftype2);    % generate process and observation noise sources for SRCDKF


%--- Some default values

pNoiseCov0 = chol(1e-4*eye(model.paramdim))';

Px0 = eye(InfDS.statedim);   % initial state covariance

%------------------- Extended Kalman Filter ------------------------------------

Xh1 = zeros(InfDS1.statedim,N);
Px1 = Px0;         % initial state covariance
pNoise1.cov(2:end,2:end) = pNoiseCov0;

fprintf('\n  EKF running... ');

[Xh1, Px1, pNoise1] = ekf(Xh1(:,1), Px1, pNoise1, oNoise1, y, [], [], InfDS1);

fprintf(' done.\n\n ');


%------------------- Square-root Unscented Kalman Filter -------------

Xh2 = zeros(InfDS1.statedim,N);

Sx2 = chol(Px0)';      % square-root filters operate on the Cholesky factor (matrix
                       % sqaure root) of the state covariance

alpha = 1e-1;
beta = 2;
kappa = 0;

InfDS2.spkfParams = [alpha beta kappa];  % UKF parameters

pNoise2.cov(2:end,2:end) = pNoiseCov0;


fprintf(' SRUKF running... ');

[Xh2, Sx2, pNoise2] = srukf(Xh2(:,1), Sx2, pNoise2, oNoise2, y, [], [], InfDS2);

fprintf(' done.\n\n ');

%---------------------------------------------------------------------------------


%--- Plot results

figure(1);
clf
subplot(211);
p1 = plot(X(1,:),'b','linewidth',2); hold on
p2 = plot(y,'g+');
p3 = plot(Xh1(1,:),'m');
p4 = plot(Xh2(1,:),'r'); hold off
legend([p1 p2 p3 p4],'clean','noisy','EKF estimate','SRUKF estimate',-1);
xlabel('time');
ylabel('state(1)');
title('DEMJE1 : LTI System Joint Estimation');

true_model_trace = cvecrep(model.params,N-1);

Wh1 = Xh1(InfDS.model.statedim+1:end,:);
Wh2 = Xh2(InfDS.model.statedim+1:end,:);

subplot(212);
p11=plot(true_model_trace(1,:),'b--','linewidth',3); hold on;
p12=plot(true_model_trace(2,:),'r--','linewidth',3);
p21=plot(Wh1(1,:),'b-.');
p22=plot(Wh1(2,:),'r-.');
p31=plot(Wh2(1,:),'b');
p32=plot(Wh2(2,:),'r'); hold off
axis([1 N -1.5 2.5]);
legend([p11 p21 p31 p12 p22 p32],'true parameter 1','EKF estimate','SRUKF estimate','true parameter 2','EKF estimate','SRUKF estimate',-1);
xlabel('time');
ylabel('parameter values');


%--- Calculate mean square estimation error

mse1 = mean((Xh1(1,:)-X(1,:)).^2);
mse2 = mean((Xh2(1,:)-X(1,:)).^2);

fprintf('\nMean-square-error (MSE) of EKF estimate   : %4.3f\n', mse1);
fprintf('\nMean-square-error (MSE) of SRUKF estimate : %4.3f\n\n', mse2);


%--- House keeping

remrelpath('../gssm');       % remove relative search path to example GSSM files from MATLABPATH
remrelpath('../data');       % remove relative search path to example data files from MATLABPATH
