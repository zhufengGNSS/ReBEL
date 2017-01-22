% DEMSE1  Demonstrate state estimation on a simple 2nd order LTI system.
%
%   This is a simple demonstration of how to use the ReBEL toolkit for state estimation on
%   a simple 2nd order LTI system.
%
%   See also
%   GSSM_LTI1
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

fprintf('\nDEMSE1 : This demonstration shows how the ReBEL toolkit is used for simple state estimation\n');
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

%ftype = input('Type of estimator [ kf, ekf, ukf, cdkf, srcdkf, or srukf ] ? ','s');

%if ~stringmatch(ftype,{'kf','ekf','ukf','cdkf','srcdkf','srukf'})
%    error('That estimator/filter type is not recognized.');
%end

%--- Use a couple of different filters...
lftype = {'kf','ekf','ukf','cdkf','srukf','srcdkf'};

for k=1:6,
    
  ftype = lftype{k};    


  %--- Setup argument data structure which serves as input to
  %--- the 'geninfds' function. This function generates the InferenceDS
  %--- data structures which are needed by all inference algorithms
  %--- in the PiLab toolkit.

  Arg.type = 'state';                                  % inference type (state estimation)
  Arg.tag = 'State estimation for GSSM_LTI1 system.';  % arbitrary ID tag
  Arg.model = model;                                   % GSSM data structure of external system

  InfDS = geninfds(Arg);                               % Create inference data structure and

  [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, ftype);       % generate process and observation noise sources

  %--- Setup runtime buffers

  Xh = zeros(InfDS.statedim,N);          % state estimation buffer
  Xh(:,1) = zeros(size(X(:,1)));         % initial estimate of state E[X(0)]
  Px = eye(InfDS.statedim);              % initial state covariance



  %--- Call inference algorithm / estimator

  switch ftype

    %------------------- Linear Kalman Filter --------------------------------------
    case 'kf'

        [Xh, Px] = kf(Xh(:,1), Px, pNoise, oNoise, y, [], [], InfDS);

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




  end


  %--- Plot results

  figure(1); subplot(3,2,k);
  p1 = plot(X(1,:)); hold on
  p2 = plot(y,'g+');
  p3 = plot(Xh(1,:),'r'); hold off;
  legend([p1 p2 p3],'clean','noisy',[ftype ' estimate']);
  xlabel('time');
  title('DEMSE1 : LTI System State Estimation');


  %--- Calculate mean square estimation error
  mse = mean((Xh(1,:)-X(1,:)).^2);
  disp([ftype ' : Mean square error (MSE) of estimate : ' num2str(mse)]);

end  
  
%--- House keeping

remrelpath('../gssm');       % remove relative search path to example GSSM files from MATLABPATH
remrelpath('../data');       % remove relative search path to example data files from MATLABPATH