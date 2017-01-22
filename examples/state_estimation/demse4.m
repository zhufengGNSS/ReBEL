% DEMSE4  Bearing Only Tracking Example
%
%   This demonstrates the use of the Sigma-Point Particle Filter on the classic
%   HARD bearing only tracking problem.
%
%   Note : This problem has not been optimised for optimal performance. It is
%          simply to demonstrate the use of the SPPF on a tracking problem.
%
%
%   See also
%   GSSM_BOT
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

fprintf('\nDEMSE4 :  Bearing Only Tracking\n\n');
fprintf('A random target trajectory is generated for each run, resulting in varying\n');
fprintf('tracking performance.\n\n');
fprintf('NOTE 1: This example has not been optimised for optimal performance.\n');
fprintf('It is used simply to demonstrate the use of the different inference\n');
fprintf('algorithms on a very hard nonlinear tracking problem.\n\n');
fprintf('NOTE 2: To best view result plots, maximize window and move legend boxes to outside frame.\n\n\n');


%--- General setup

addrelpath('../gssm');         % add relative search path to example GSSM files to MATLABPATH
addrelpath('../data');         % add relative search path to example data files to MATLABPATH

%--- Initialise GSSM

model = gssm_bot('init');

%--- Generate inference data structure

Arg.model = model;                                   % embed GSSM
Arg.type = 'state';                                  % estimation type
Arg.tag = 'State estimation for bearings-only tracking problem';  % info tag (not required)

InfDS = geninfds(Arg);                               % call generate function

%--- Generate estimation process and observation noise sources
%ftype = input('Inference algorithm  [ srcdkf / pf / sppf / gmsppf ] : ','s');  %  set type of inference algorithm (estimator) to use :


%--- Generate some data : Initial target state generated according to Gordon, Salmond & Ewing - 1995

N = 25;                                              % max. time k=1..N

V = model.pNoise.sample( model.pNoise, N);     % generate process noise
W = model.oNoise.sample( model.oNoise, N);     % generate observation noise

X = zeros(InfDS.statedim, N);                        % system state buffer
y = zeros(InfDS.obsdim,N);                           % system observations buffer

bearing_0      = -pi+rand(1)*2*pi;
bearing_rate_0 = 0.1*randn(1);
range_0        = 0.1*randn(1)+1;
range_rate_0   = 0.01*randn(1)-0.1;

X(:,1) = [range_0*cos(bearing_0);                       % initial target location in 2D-cartesian space
         (range_0 + range_rate_0)*cos(addangle(bearing_0,bearing_rate_0)) - range_0*cos(bearing_0);
         range_0*sin(bearing_0);
         (range_0 + range_rate_0)*sin(addangle(bearing_0,bearing_rate_0)) - range_0*sin(bearing_0)];

y(:,1) = model.hfun( model, X(:,1), W(:,1), []);  % initial observation

for k=2:N,
    X(:,k) = model.ffun( model, X(:,k-1), V(:,k-1), []);
    y(:,k) = model.hfun( model, X(:,k), W(:,k), []);
end

true_range   = sqrt(X(1,:).^2 + X(3,:).^2);             % calculate range ground truth trajectory
true_bearing = atan2(X(3,:), X(1,:));                   % calculate bearing ground truth trajectory

%--- Use a couple of different filters...
lftype = {'srcdkf','pf','sppf','gmsppf'};

figure(1); clf;

for k=1:4,
    
  ftype = lftype{k};    

  %--- Setup estimation buffers

  Xh = zeros(InfDS.statedim, N);
  Sx = eye(InfDS.statedim);

  range_error   = zeros(1,N);
  bearing_error = zeros(1,N);
  pos_error     = zeros(1,N);

  %--- Determine initial uncertainty in vehicle position

  Nstat = 10000;

  Wstat = model.oNoise.sample( model.oNoise, Nstat);

  bearing_stat      = bearing_0 + sqrt(model.oNoise.cov(1,1))*randn(1,Nstat);
  bearing_rate_stat = 0.1*randn(1,Nstat);
  range_stat        = 0.1*randn(1,Nstat)+1;
  range_rate_stat   = 0.01*randn(1,Nstat)-0.1;

  Xstat = [range_stat.*cos(bearing_stat);
         (range_stat + range_rate_stat).*cos(addangle(bearing_stat,bearing_rate_stat)) - range_stat.*cos(bearing_stat);
         range_stat.*sin(bearing_stat);
         (range_stat + range_rate_stat).*sin(addangle(bearing_stat,bearing_rate_stat)) - range_stat.*sin(bearing_stat)];

  Mu0 = mean(Xstat,2);
  P0  = cov(Xstat');

  Xh(:,1) = Mu0;                  % initial state distribution : mean
  Sx = chol(P0)';                 % initial state distribution : covariance Cholesky factor


  %--- Display target trajectory detail

  figure(1); subplot(3,4,k);
  p1=plot(X(1,:),X(3,:),'-*'); hold on;
  p2=plot(X(1,1),X(3,1),'c*');
  p3=plot(X(1,end),X(3,end),'m*');
  p4=plot(0,0,'kx','linewidth',2); hold off;
  %legend([p1 p2 p3 p4],'target trajectory','position : k=0',['position : k=' num2str(N)],'observer position',0);
  xlabel('x');
  ylabel('y');
  title(['Target Trajectory - ' ftype]);
  axis tight
  vmin1 = axis;
  vmax1 = axis;
  
  subplot(3,4,k+4);
  p11=plot(1:N,true_range,'b-o');
  xlabel('k');
  ylabel('range');
  title(['Range Profile - ' ftype]);
  %legend([p11],'true');
  axis tight
  vmin2 = axis;
  vmax2 = axis;
  subplot(3,4,k+8);
  p13=plot(1:N,true_bearing,'b-o'); hold on;
  p14=plot(1:N,y(1,:),'k+'); hold off;
  %legend([p13 p14],'true bearing','measured bearing',0);
  xlabel('time : k');
  ylabel('bearing : radians');
  title(['Bearing Profile - ' ftype]);
  axis tight
  vmin3 = axis;
  vmax3 = axis;
  
  drawnow


%-------------------------------------------------------

switch ftype
case {'pf','gspf','gmsppf'}
  numParticles = 1000;                        % number of particles
otherwise
  numParticles = 200;
end

bearing_stat      = bearing_0+sqrt(model.oNoise.cov(1,1))*randn(1,numParticles);
bearing_rate_stat = 0.1*randn(1,numParticles);
range_stat        = 0.1*randn(1,numParticles)+1;
range_rate_stat   = 0.01*randn(1,numParticles)-0.1;

initialParticles = [range_stat.*cos(bearing_stat);
                    (range_stat + range_rate_stat).*cos(addangle(bearing_stat,bearing_rate_stat)) - range_stat.*cos(bearing_stat);
                    range_stat.*sin(bearing_stat);
                    (range_stat + range_rate_stat).*sin(addangle(bearing_stat,bearing_rate_stat)) - range_stat.*sin(bearing_stat)];

initialParticles = Sx*randn(InfDS.statedim,numParticles) + cvecrep(Mu0,numParticles);

initialParticlesCov = repmat(Sx,[1 1 numParticles]);  % particle covariances


%=================================================================================================================
%--- Run estimator on observed data (noisy bearing readings)

disp([ftype ' : Estimating trajectory...']);

  switch ftype

  %---------------------------------------------------------------------------------------------------------
  case 'pf'

      [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, ftype);     % call system noise sources generation function

      ParticleFiltDS.N = numParticles;
      ParticleFiltDS.particles = initialParticles;
      ParticleFiltDS.weights = (1/numParticles)*ones(1,numParticles);
      InfDS.resampleThreshold = 1;    % set resample threshold
      InfDS.estimateType = 'mean';    % estimate type for Xh

      [Xh, ParticleFiltDS] = pf(ParticleFiltDS, pNoise, oNoise, y, [], [], InfDS);

  %---------------------------------------------------------------------------------------------------------
  case 'gmsppf'

      [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, ftype);     % call system noise sources generation function

      ParticleFiltDS.N = numParticles;      % number of particles
      ParticleFiltDS.stateGMM = gmmfit(initialParticles, 5, [0.001 10], 'sqrt');  % fit a 5 component GMM to initial state distribution
      InfDS.estimateType = 'mean';    % estimate type for Xh
      InfDS.spkfType = 'srcdkf';      % Type of SPKF to use inside SPPF (note that ParticleFiltDS.particlesCov should comply)
      InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)

      [Xh, ParticleFiltDS] = gmsppf(ParticleFiltDS, pNoise, oNoise, y, [], [], InfDS);


  %---------------------------------------------------------------------------------------------------------
  case 'sppf'

      [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, ftype);     % call system noise sources generation function

      InfDS.spkfType = 'srcdkf';      % Type of SPKF to use inside SPPF (note that ParticleFiltDS.particlesCov should comply)
      InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)
      InfDS.resampleThreshold = 1;    % set resample threshold
      InfDS.estimateType = 'mean';    % estimate type for Xh

      [pNoiseGAUS, oNoiseGAUS, foo] = gensysnoiseds(InfDS, InfDS.spkfType); % generate Gaussian system noise sources for internal SPKFs

      % build ParticleFilter data structure

      ParticleFiltDS.N = numParticles;              % number of particles
      ParticleFiltDS.particles = initialParticles;  % initialize particle means
      ParticleFiltDS.particlesCov = initialParticlesCov;  % initialize article covariances
      ParticleFiltDS.pNoise = pNoiseGAUS;      % embed SPKF noise sources
      ParticleFiltDS.oNoise = oNoiseGAUS;      %   "   "       "    "
      ParticleFiltDS.weights = cvecrep(1/numParticles,numParticles); % initialize particle weights

      [Xh, ParticleFiltDS] = sppf(ParticleFiltDS, pNoise, oNoise, y, [], [], InfDS);


  %---------------------------------------------------------------------------------------------------------
  case 'srcdkf'

      [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, ftype);     % call system noise sources generation function

      InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)

      [Xh, Sx] = srcdkf(Xh(:,1), Sx, pNoise, oNoise, y, [], [], InfDS);



  %---------------------------------------------------------------------------------------------------------
  case 'srukf'

      [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, ftype);     % call system noise sources generation function


      InfDS.spkfParams  = [1 2 0];    % scale factor (CDKF parameter h)


      [Xh, Sx] = srukf(Xh(:,1), Sx, pNoise, oNoise, y, [], [], InfDS);



  %---------------------------------------------------------------------------------------------------------
  otherwise

      error([' Unknown inference algorithm type ''' ftype '''']);


end


%=================================================================================================================


%--- Calculate errors

range_estimate = sqrt(Xh(1,:).^2 + Xh(3,:).^2);
bearing_estimate = atan2(Xh(3,:), Xh(1,:));

range_error   =  range_estimate - true_range;
bearing_error =  bearing_estimate - true_bearing;
pos_error     =  sqrt((Xh([1;3],:)-X([1;3],:)).^2);


%--- Display results

figure(1); subplot(3,4,k); hold on;
p5=plot(Xh(1,:),Xh(3,:),'r-o');
plot(Xh(1,1),Xh(3,1),'c*');
plot(Xh(1,end),Xh(3,end),'m*');
if (k==1), legend([p1 p2 p3 p4 p5],'trajectory','position: k=0',['position: k=' num2str(N)],'observer','estimate',0); end
xlabel('x');
ylabel('y');
title(['Target Trajectory - ' ftype]);
axis tight
vmin1 = min([vmin1; axis]);
vmax1 = max([vmax1; axis]);
hold off;

figure(1); 
subplot(3,4,k+4); hold on;
p12=plot(1:N,range_estimate,'r-');
xlabel('k');
ylabel('range');
title(['Range Profile - ' ftype]);
%legend([p11 p12],'true','inferred');
axis tight
vmin2 = min([vmin2; axis]);
vmax2 = max([vmax2; axis]);

hold off;
subplot(3,4,k+8); hold on
p15=plot(1:N,bearing_estimate,'r-');
xlabel('t');
ylabel('bearing');
title(['Bearing Profile - ' ftype])
axis tight
vmin3 = min([vmin3; axis]);
vmax3 = max([vmax3; axis]);
if (k==1), legend([p13 p14 p15],'true','measured','inferred',0); end
hold off;

end

for k=1:4,
    figure(1);
    subplot(3,4,k); axis([vmin1(1) vmax1(2) vmin1(3) vmax1(4)]);
    subplot(3,4,k+4); axis([vmin2(1) vmax2(2) vmin2(3) vmax2(4)]);
    subplot(3,4,k+8); axis([vmin3(1) vmax3(2) vmin3(3) vmax3(4)]);
end

%--- House keeping

remrelpath('../gssm');       % remove relative search path to example GSSM files from MATLABPATH
remrelpath('../data');       % remove relative search path to example data files from MATLABPATH
