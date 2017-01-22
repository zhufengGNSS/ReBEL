% DEMSE2  Demonstrate state estimation on a simple scalar nonlinear (time variant) problem
%
%   See also
%   GSSM_N1
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

fprintf('\nDEMSE2 : This demonstration shows how the ReBEL toolkit is used for simple state estimation\n');
fprintf('         on a scalar nonlinear problem.\n\n');


%--- General setup

addrelpath('../gssm');         % add relative search path to example GSSM files to MATLABPATH
addrelpath('../data');         % add relative search path to example data files to MATLABPATH
addrelpath('../../netlab');    % Some of the algorithms requires Netlab functions

%--- Ask the user which inference algorithm to use
%ftype = input('Type of estimator [ ekf, ukf, cdkf, srukf, srcdkf, pf, gspf, gmsppf or sppf ] ? ','s');
%if ~stringmatch(ftype,{'ekf','cdkf','ukf','srukf','srcdkf','pf','gspf','sppf','gmsppf'})
%    error('That estimator/filter type is not recognized.');
%end

%--- Compare these algorithms...
lftype={'ekf','cdkf','ukf','srukf','srcdkf','pf','gspf','sppf','gmsppf'};

number_of_runs = input('Number of independent runs ? ');

mean_RMSE = zeros(1,9);   % buffer for MC results for each algorithm
var_RMSE  = zeros(1,9);   %     "                              "

for jj=1:9,

  ftype=lftype{jj};

  disp(['[' ftype ']']);

  %--- Initialise GSSM model from external system description script.

  model = gssm_n1('init');

  Arg.type = 'state';                                  % inference type (state estimation)
  Arg.tag = 'State estimation for GSSM_N1 system.';    % arbitrary ID tag
  Arg.model = model;                                   % GSSM data structure of external system
  Arg.algorithm = ftype;                               % set inference algorithm to be used

  InfDS = geninfds(Arg);                               % Create inference data structure and
  [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS,ftype);       % generate process and observation noise sources


  %--- Loop over number of independent runs

  for k=1:number_of_runs

    randn('state',sum(100*clock));          % stir the pot... shuffle the deck :-)
    rand('state',sum(100*clock));

    %--- Generate some data

    N  = 60;                                                % number of datapoints
    X  = zeros(model.statedim,N);                           % state data buffer
    y  = zeros(model.obsdim,N);                             % observation data buffer

    pnoise = model.pNoise.sample( model.pNoise, N);   % generate process noise
    onoise = model.oNoise.sample( model.oNoise, N);   % generate observation noise

    X(1) = 1;                                               % initial state
    y(1) = model.hfun( model, X(1), onoise(1), 1);    % observation of initial state
    for j=2:N,
      X(j) = model.ffun( model, X(:,j-1), pnoise(j-1), j-1);
      y(j) = model.hfun( model, X(:,j), onoise(j), j);
    end

    U1 = [0:N-1];
    U2 = [1:N];

    %--- Setup runtime buffers

    Xh = zeros(InfDS.statedim,N);          % state estimation buffer
    Xh(:,1) = 1;                           % initial estimate of state E[X(0)]
    Px = 3/4*eye(InfDS.statedim);          % initial state covariance

    %--- Call inference algorithm / estimator

    switch ftype


      %------------------- Extended Kalman Filter ------------------------------------
      case 'ekf'

        [Xh, Px] = ekf(Xh(:,1), Px, pNoise, oNoise, y, U1, U2, InfDS);


      %------------------- Unscented Kalman Filter -----------------------------------
      case 'ukf'

        alpha = 1;         % scale factor (UKF parameter)
        beta  = 2;         % optimal setting for Gaussian priors (UKF parameter)
        kappa = 0;         % optimal for state dimension=2 (UKF parameter)

        InfDS.spkfParams = [alpha beta kappa];

        [Xh, Px] = ukf(Xh(:,1), Px, pNoise, oNoise, y, U1, U2, InfDS);


      %------------------- Central Difference Kalman Filter ---------------------------
      case 'cdkf'

        InfDS.spkfParams = sqrt(3);    % scale factor (CDKF parameter h)

        [Xh, Px] = cdkf(Xh(:,1), Px, pNoise, oNoise, y, U1, U2, InfDS);


      %------------------- Square Root Unscented Kalman Filter ------------------------
      case 'srukf'

        alpha = 1;         % scale factor (UKF parameter)
        beta  = 2;         % optimal setting for Gaussian priors (UKF parameter)
        kappa = 0;         % optimal for state dimension=2 (UKF parameter)

        Sx = chol(Px)';

        InfDS.spkfParams = [alpha beta kappa];

        [Xh, Sx] = srukf(Xh(:,1), Sx, pNoise, oNoise, y, U1, U2, InfDS);


      %------------------- Square Root Central Difference Kalman Filter ---------------
      case 'srcdkf'

        InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)

        Sx = chol(Px)';

        [Xh, Sx] = srcdkf(Xh(:,1), Sx, pNoise, oNoise, y, U1, U2, InfDS);


      %------------------- Generic Particle Filter (a.k.a Bootstrap-filter of CONDENSATION -----------
      case 'pf'

        M = 200;                             % number of particles
        ParticleFiltDS.N = M;
        ParticleFiltDS.particles = randn(InfDS.statedim,M)+cvecrep(Xh(:,1),M);  % initialize particles
        ParticleFiltDS.weights = cvecrep(1/M,M); % initialize weights

        InfDS.resampleThreshold = 0.5;    % set resample threshold
        InfDS.estimateType = 'mean';      % estimate type for Xh

        [Xh, ParticleFiltDS] = pf(ParticleFiltDS, pNoise, oNoise, y, U1, U2, InfDS);

      %------------------- Gaussian-Sum Particle Filter ---------------------------------------------
      case 'gspf'

        M = 200;                             % number of particles
        ParticleFiltDS.N = M;

        initialParticles = randn(InfDS.statedim,M)+cvecrep(Xh(:,1),M);  % initialize particles

        ParticleFiltDS.stateGMM = gmmfit(initialParticles, 2, [0.001 10], 'sqrt');  % fit a 3 component GMM to initial state distribution

        InfDS.estimateType = 'mean';      % estimate type for Xh
        InfDS.threshold = 0.001;

        Arg.type='gmm';
        Arg.cov_type='sqrt';
        Arg.dim=model.Vdim;
        Arg.M = 2;
        Arg.mu = cvecrep(model.pNoise.mu,Arg.M);
        Arg.cov = zeros(Arg.dim,Arg.dim,Arg.M);
        Arg.cov(:,:,1) = 2*model.pNoise.cov(:,:,1);
        Arg.cov(:,:,2) = 0.5*model.pNoise.cov(:,:,1);
        Arg.weights = [0.5 0.5];
        pNoise = gennoiseds(Arg);

        [Xh, ParticleFiltDS] = gspf(ParticleFiltDS, pNoise, oNoise, y, U1, U2, InfDS);

      %------------------- Sigma-Point Bayes Filter ---------------------------------------------
      case 'gmsppf'

        M = 200;
        ParticleFiltDS.N = M;            % number of particles

        initialParticles = randn(InfDS.statedim,M)+cvecrep(Xh(:,1),M);  % initialize particles

        tempCov = zeros(1,1,2); tempCov(:,:,1) = sqrt(2); tempCov(:,:,2)=1;

        ParticleFiltDS.stateGMM = gmmfit(initialParticles, 3, [0.001 10], 'sqrt');  % fit a 3 component GMM to initial state distribution

        InfDS.estimateType = 'mean';    % estimate type for Xh

        InfDS.spkfType = 'srcdkf';      % Type of SPKF to use inside SPPF (note that ParticleFiltDS.particlesCov should comply)
        InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)

        Arg.type='gmm';
        Arg.cov_type='sqrt';
        Arg.dim=model.Vdim;
        Arg.M = 2;
        Arg.mu = cvecrep(model.pNoise.mu,Arg.M);
        Arg.cov = zeros(Arg.dim,Arg.dim,Arg.M);
        Arg.cov(:,:,1) = 2*model.pNoise.cov(:,:,1);
        Arg.cov(:,:,2) = 0.5*model.pNoise.cov(:,:,1);
        Arg.weights = [0.5 0.5];
        pNoise = gennoiseds(Arg);

        [Xh, ParticleFiltDS] = gmsppf(ParticleFiltDS, pNoise, oNoise, y, U1, U2, InfDS);



      %------------------- Sigma-Point Particle Filter -----------------------------------------------
      case 'sppf'

        M = 200;                             % number of particles
        ParticleFiltDS.N = M;
        ParticleFiltDS.particles  = cvecrep(Xh(:,1),M);  % initialize particle means
        ParticleFiltDS.particlesCov = repmat(eye(InfDS.statedim),[1 1 M]);       % particle covariances

        pNoiseGAUS.cov = sqrt(2*3/4);
        oNoiseGAUS.cov = sqrt(1e-1);

        [pNoiseGAUS, oNoiseGAUS, foo] = gensysnoiseds(InfDS,'srukf');

        ParticleFiltDS.pNoise = pNoiseGAUS;
        ParticleFiltDS.oNoise = oNoiseGAUS;
        ParticleFiltDS.weights = cvecrep(1/M,M); % initialize weights

        InfDS.spkfType = 'srukf';         % Type of SPKF to use (note that ParticleFiltDS.particlesP should comply)
        %InfDS.spkfParams = [sqrt(3)];
        InfDS.spkfParams = [1 0 2];
        InfDS.resampleThreshold = 1;    % set resample threshold
        InfDS.estimateType = 'mean';      % estimate type for Xh

        [Xh, ParticleFiltDS] = sppf(ParticleFiltDS, pNoise, oNoise, y, U1, U2, InfDS);

    end

    %--- Plot results

    figure(1); clf;
    p1 = plot(X(1,:)); hold on
    p2 = plot(y,'g+');
    p3 = plot(Xh(1,:),'r'); hold off;
    legend([p1 p2 p3],'clean','noisy',[ftype ' estimate']);
    xlabel('time');
    title('DEMSE2 : Nonlinear Time Variant State Estimation (non Gaussian noise)');

    drawnow

    %--- Calculate mean square estimation error

    rmse(k) = sqrt(mean((Xh(1,2:end)-X(1,2:end)).^2));
    fprintf('%d:%d  Root-mean-square-error (RMSE) of estimate : %4.3f\n', k, number_of_runs, rmse(k));

  end

  disp(' ');
  mean_RMSE(jj) = mean(rmse);
  var_RMSE(jj)  = var(rmse);

end


%--- Summary of results
disp(' ');
disp(' ');
for jj=1:9,
   disp([lftype{jj} ' - RMSE (mean) : ' num2str(mean_RMSE(jj)) '   RMSE (var) : ' num2str(var_RMSE(jj))]);
end

%--- House keeping

remrelpath('../gssm');       % remove relative search path to example GSSM files from MATLABPATH
remrelpath('../data');       % remove relative search path to example data files from MATLABPATH
