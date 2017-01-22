% GSSM_LTI1  Generalized state space model for simple LTI system
%
% The model is a simple 2nd order LTI system with Gaussian observation noise
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

  model.type = 'gssm';                       % object type = generalized state space model
  model.tag  = 'GSSM_LTI1';                  % ID tag

  model.ffun       = @ffun;                  % file handle to FFUN
  model.hfun       = @hfun;                  % file handle to HFUN
  model.prior      = @prior;                 % file handle to PRIOR
  model.likelihood = @likelihood;            % file handle to LIKELIHOOD
  model.linearize  = @linearize;             % file handle to LINEARIZE
  model.setparams  = @setparams;             % file handle to SETPARAMS

  model.statedim   = 2;                      %   state dimension
  model.obsdim     = 1;                      %   observation dimension
  model.paramdim   = 2;                      %   parameter dimension
  model.U1dim      = 0;                      %   exogenous control input 1 dimension
  model.U2dim      = 0;                      %   exogenous control input 2 dimension
  model.Vdim       = 1;                      %   process noise dimension
  model.Ndim       = 1;                      %   observation noise dimension

  Arg.type = 'gaussian';
  Arg.cov_type = 'full';
  Arg.dim = model.Vdim;
  Arg.mu = 0;
  Arg.cov  = 0.001;
  model.pNoise     = gennoiseds(Arg);   % process noise : zero mean white Gaussian noise , cov=0.001

  Arg.type = 'gaussian';
  Arg.cov_type = 'full';
  Arg.dim = model.Ndim;
  Arg.mu = 0;
  Arg.cov  = 0.3;
  model.oNoise     = gennoiseds(Arg);     % observation noise : zero mean white Gaussian noise, cov=0.2

  model.params     = zeros(model.paramdim,1);
  model.A  = [model.params(1) model.params(2); 1 0];
  model.B  = [];
  model.C  = [1 0];
  model.D  = [];
  model.G  = [1 0]';
  model.H  = [1];

  model = setparams(model,[1.9223 -0.9604]);   % 2nd order under-damped LTI system


%===============================================================================================
function model = setparams(model, params, index_vector)

  if (nargin==2)
    model.params = params(:);
  elseif (nargin==3)
    model.params(index_vector) = params(:);
  else
    error('[ setparams ] Incorrect number of input arguments.');
  end

  model.A(1,:)  = model.params';

%===============================================================================================
function new_state = ffun(model, state, V, U1)

  new_state      = model.A * state;

  if ~isempty(V)
      new_state(1,:) = new_state(1,:) + V(1,:);
  end


%===============================================================================================
function tranprior = prior(model, nextstate, state, U1, pNoiseDS)

  X = nextstate - ffun(model, state, [], U1);

  tranprior = pNoiseDS.likelihood( pNoiseDS, X(1,:));

%===============================================================================================
function observ = hfun(model, state, N, U2)

  observ = state(1,:);

  if ~isempty(N)
    observ = state(1,:) + N(1,:);
  end

%===============================================================================================
function llh = likelihood(model, obs, state, U2, oNoiseDS)

  X = obs - hfun(model, state, [], U2);

  llh = oNoiseDS.likelihood( oNoiseDS, X);


%===============================================================================================
function out = linearize(model, state, V, N, U1, U2, term, index_vector)

  if (nargin<7)
    error('[ linearize ] Not enough input arguments!');
  end

  %--------------------------------------------------------------------------------------
  switch (term)

    case 'A'
      %%%========================================================
      %%%             Calculate A = df/dstate
      %%%========================================================
      out = model.A;

    case 'B'
      %%%========================================================
      %%%             Calculate B = df/dU1
      %%%========================================================
      out = model.B;

    case 'C'
      %%%========================================================
      %%%             Calculate C = dh/dx
      %%%========================================================
      out = model.C;

    case 'D'
      %%%========================================================
      %%%             Calculate D = dh/dU2
      %%%========================================================
      out = model.D;

    case 'G'
      %%%========================================================
      %%%             Calculate G = df/dv
      %%%========================================================
      out = model.G;

    case 'H'
      %%%========================================================
      %%%             Calculate H = dh/dn
      %%%========================================================
      out = model.H;

    case 'JFW'
      %%%========================================================
      %%%             Calculate  = dffun/dparameters
      %%%========================================================
      out = [state(1) state(2); 0 0];


    case 'JHW'
      %%%========================================================
      %%%             Calculate  = dhfun/dparameters
      %%%========================================================
      out = zeros(model.obsdim,model.paramdim);

    otherwise
      error('[ linearize ] Invalid model term requested!');

  end

  if (nargin==8), out = out(:,index_vector); end

  %--------------------------------------------------------------------------------------