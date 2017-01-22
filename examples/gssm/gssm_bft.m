% GSSM_BFT  General state space model for Bearings and Frequency Tracking of a randomly maneuvering
%           target relative to a moving observer (submarine).
%
%   The following state space model is used :
%
%     X(k) = |1 4 0 0 0| X(k-1) + | 8  0  0| V(k-1)
%            |0 1 0 0 0|          | 4  0  0|
%            |0 0 1 4 0|          | 0  8  0|
%            |0 0 0 1 0|          | 0  4  0|
%            |0 0 0 0 1|          | 0  0  1|
%    (the measures are given by the sonar every 4 seconds)
%
%   Where the state vector is defined as
%            - the 2D position and velocity vector of the target (relative to a fixed external reference frame),
%            - the pure tone frequency emitted by the target (very stable).
%
%     X(k) = |x1(k)| = |      x-position at time k     |
%            |x2(k)|   |      x-velocity at time k     |
%            |x3(k)|   |      y-position at time k     |
%            |x4(k)|   |      y-velocity at time k     |
%            |x5(k)|   | pure tone frequency at time k |
%
%   And the observations at time k, O(k) are :
%           - the bearing angle (in radians) from the moving observer (submarine) towards the target,
%           - the doppler-shifted frequency tone tracked by the observer (submarine).
%
%     O(k) = |                bearing = arctan((x3(k)-sub3(k))/(x1(k)-sub1(k)))                            |  +  | v1(k) |
%            |   frequency = x5(k)*(1+1/1500*((x2(k)-sub2(k))*cos(bearing)+(x4(k)-sub4(k))*sin(bearing)))  |     | v2(k) |
%
%   c=1500 m/s (sound speed in water)
%
%   The submarine state is known at each time precisely and described by the following vector :
%     sub(k) = |sub1(k)| = |  x-position at time k     |
%              |sub2(k)|   |  x-velocity at time k     |
%              |sub3(k)|   |  y-position at time k     |
%              |sub4(k)|   |  y-velocity at time k     |
%              |sub5(k)|   |  frequency tone at time k | (not used here)
%
%   The state dynamics are driven by a 2 dimensional white Gaussian noise source and the
%   observations are corrupted by additive scalar white Gaussian noise.
%
%   See :  Gordon, Salmond & Ewing, "Bayesian State Estimation for Tracking and Guidance Using
%   the Bootstrap Filter", Journal of Guidance, Control and Dynamics, 1995.
%
%   This example is courtesy of Alain Bonnot.
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

function [varargout] = model_interface(func, varargin)

  switch func

    %--- Initialize GSSM data structure --------------------------------------------------------
    case 'init'
      model = init(varargin);
        error(consistent(model,'gssm'));               % check consistentency of initialized model
      varargout{1} = model;

    %--------------------------------------------------------------------------------------------
    otherwise

      error(['Function ''' func ''' not supported.']);

  end


%===============================================================================================
function model = init(init_args)

  model.type = 'gssm';                              % object type = generalized state space model
  model.tag  = 'GSSM_Bearings_Frequency_Tracking';  % ID tag

  model.statedim   = 5;                      %   state dimension
  model.obsdim     = 2;                      %   observation dimension
  model.paramdim   = 12;                     %   parameter dimension
                                             %   parameter estimation will be done)
  model.U1dim      = 0;                      %   exogenous control input 1 dimension
  model.U2dim      = 5;                      %   exogenous control input 2 dimension
  model.Vdim       = 3;                      %   process noise dimension
  model.Ndim       = 2;                      %   observation noise dimension

  model.ffun      = @ffun;                   % file handle to FFUN
  model.hfun      = @hfun;                   % file handle to HFUN
  model.prior     = @prior;
  model.likelihood = @likelihood;            % file handle to LIKELIHOOD
  model.innovation = @innovation;            % file handle to INNOVATION
  model.setparams  = @setparams;             % file handle to SETPARAMS

  model.obsAngleCompIdxVec = [1];            % indicate that the first (and only component) of the observation
                                             % vector is an angle measured in radians. This is needed so that the
                                             % SPKF based algorithms can correctly deal with the angular discontinuity
                                             % at +- pi radians.


  Arg.type = 'gaussian';
  Arg.dim = model.Vdim;
  Arg.mu  = zeros(Arg.dim,1);
  Arg.cov_type = 'full';
  Arg.cov   = [((1e-3)^2)*eye(Arg.dim-1) zeros(Arg.dim-1,1);zeros(1,Arg.dim-1) (1e-4)^2];
  model.pNoise = gennoiseds(Arg);            % process noise : zero mean white Gaussian noise, cov = (1e-3)^2(dynamics) and (1e-4)^2 (tone frequency, very stable)

  Arg.type = 'gaussian';
  Arg.dim = model.Ndim;
  Arg.mu = zeros(Arg.dim,1);
  Arg.cov_type ='full';
  Arg.cov  = [0.0175^2 0;0 0.06^2];
  model.oNoise = gennoiseds(Arg);            % observation noise : zero mean white Gaussian noise, cov=0.0175^2 (bearings=1° error) and 0.06^2 (frequencies)

  model.params = zeros(model.paramdim,1);
  model.A = zeros(model.statedim, model.statedim);
  model.G = zeros(model.statedim, model.Vdim);

  model = setparams(model,[1 4 1 1 4 1 1 8 4 8 4 1]');


%===============================================================================================
%-- Unpack and update model internal parameters from parameter vector, 'params'

function model = setparams(model, params, index_vector)

  if (nargin==2)
    model.params = params(:);
  elseif (nargin==3)
    model.params(index_vector) = params(:);
  else
    error('[ setparams ] Incorrect number of input arguments.');
  end

  model.A([1 6 7 13 18 19 25]) = params(1:7);
  model.G([1 2 8 9 15]) = params(8:12);

  G = model.G;

  model.convFact1 = (G'*G)\(G');    % conversion matrix needed to calculate state transition prior


%===============================================================================================
%-- State transition function (vehicle dynamic model)

function new_state = ffun(model, state, V, U1)

  if isempty(V)
      new_state = model.A*state;
  else
      new_state = model.A*state + model.G*V;
  end


%===============================================================================================
%-- State observation function

function observ = hfun(model, state, N, U2)

  observ_ = zeros(2,size(state,2));

  observ_(1,:) = atan2(state(3,:)-U2(3,:),state(1,:)-U2(1,:));
  observ_(2,:) = state(5,:).*(1+1/1500*((U2(2,:)-state(2,:)).*cos(observ_(1,:))+(U2(4,:)-state(4,:)).*sin(observ_(1,:))));

%  for m=1:size(state,2),
%    observ(1,m) = atan2(state(3,m)-U2(3,m),state(1,m)-U2(1,m));
%    observ(2,m) = state(5,m)*(1+1/1500*((U2(2,m)-state(2,m))*cos(observ(1,m))+(U2(4,m)-state(4,m))*sin(observ(1,m))));
%  end


  % Now add the measurement noise... taking care with the discontinueties at +-pi radians
  if isempty(N),
    observ = observ_;
  else
      observ = observ_ + N;
      observ(1,:) = addangle(observ_(1,:), N(1,:));
  end


 %===============================================================================================
function tranprior = prior(model, nextstate, state, U1, pNoiseDS)

  V = model.convFact1 * (nextstate - model.A*state);

  tranprior = pNoiseDS.likelihood( pNoiseDS, V);

%===============================================================================================
function llh = likelihood(model, obs, state, U2, oNoiseDS)

  observ =  hfun(model, state, [], U2);

  N = obs - observ;
  N(1,:) = subangle(obs(1,:), observ(1,:));

  % Calculate log likelihood
  llh = oNoiseDS.likelihood( oNoiseDS, N);



%===============================================================================================
function innov = innovation(model, obs, observ)

  innov = obs - observ;

  % deal with the discontinueties at +-pi radians
  innov(1,:) = subangle(obs(1,:),observ(1,:));