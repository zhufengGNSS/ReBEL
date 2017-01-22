% GSSM_BOT  General state space model for Bearings-Only Tracking of a randomly maneuvering
%           target relative to a stationary observer.
%
%   The following state space model is used :
%
%     X(k) = |1 1 0 0| X(k-1) + |0.5  0 | V(k-1)
%            |0 1 0 0|          | 1   0 |
%            |0 0 1 1|          | 0  0.5|
%            |0 0 0 1|          | 0   1 |
%
%     O(k) = arctan(x3(k)/x1(k)) + n(k)
%
%   Where the state vector is defined as the 2D position and velocity vector of the target,
%   relative to a fixed external reference frame, i.e.
%
%     X(k) = |x1(k)| = |x-position at time k|
%            |x2(k)|   |x-velocity at time k|
%            |x3(k)|   |y-position at time k|
%            |x4(k)|   |y-velocity at time k|
%
%   and the observation at time k, O(k) is the bearing angle (in radians) from the fixed
%   observer towards the target.
%
%   The state dynamics are driven by a 2 dimensional white Gaussian noise source and the
%   observations are corrupted by additive scalar white Gaussian noise.
%
%   See :  Gordon, Salmond & Ewing, "Bayesian State Estimation for Tracking and Guidance Using
%   the Bootstrap Filter", Journal of Guidance, Control and Dynamics, 1995.
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

  model.type = 'gssm';                         % object type = generalized state space model
  model.tag  = 'GSSM_Bearings_Only_Tracking';  % ID tag

  model.statedim   = 4;                      %   state dimension
  model.obsdim     = 1;                      %   observation dimension
  model.paramdim   = 10;                     %   parameter dimension
                                             %   parameter estimation will be done)
  model.U1dim      = 0;                      %   exogenous control input 1 dimension
  model.U2dim      = 0;                      %   exogenous control input 2 dimension
  model.Vdim       = 2;                      %   process noise dimension
  model.Ndim       = 1;                      %   observation noise dimension

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
  Arg.cov_type = 'sqrt';
  Arg.dim = model.Vdim;
  Arg.mu  = zeros(Arg.dim,1);
  Arg.cov   = 0.01*eye(Arg.dim);
  model.pNoise = gennoiseds(Arg);            % process noise : zero mean white Gaussian noise, cov = 0.001^2

  Arg.type = 'gaussian';
  Arg.cov_type = 'sqrt';
  Arg.dim = model.Ndim;
  Arg.mu = 0;
  Arg.cov  = 0.01;
  model.oNoise = gennoiseds(Arg);            % observation noise : zero mean white Gaussian noise, cov=0.01^2

  model.params = zeros(model.paramdim,1);
  model.A = zeros(model.statedim, model.statedim);
  model.G = zeros(model.statedim, model.Vdim);

  model = setparams(model,[1 1 1 1 1 1 0.5 1 0.5 1]');


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

  model.A([1 5 6 11 15 16]) = params(1:6);
  model.G([1 2 7 8]) = params(7:10);

  G = model.G;

  model.convFact1 = (G'*G)\G';    % conversion matrix needed to calculate state transition prior


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

  observ = atan2(state(3,:),state(1,:));

  % Now add the process noise... taking care with the discontinueties at +-pi radians
  if ~isempty(N),
      observ = observ + N;
      idx = find(abs(observ) > pi);
      temp = rem(observ(idx),2*pi);
      observ(idx) = temp - sign(temp).*(2*pi);
  end


%===============================================================================================
function tranprior = prior(model, nextstate, state, U1, pNoiseDS)

  V = model.convFact1 * (nextstate - model.A*state);

  tranprior = pNoiseDS.likelihood( pNoiseDS, V);


%===============================================================================================
function llh = likelihood(model, obs, state, U2, oNoiseDS)

  observ = hfun(model, state, [], U2);

  N = subangle(obs, observ);

  llh = oNoiseDS.likelihood( oNoiseDS, N);



%===============================================================================================
function innov = innovation(model, obs, observ)

  innov = subangle(obs,observ);

