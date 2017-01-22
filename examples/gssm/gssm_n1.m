% GSSM_N1  Generalized state space model for simple nonlinear system
%
% The model is a simple scalar nonlinear system with Gamma process and Gaussian observation noise.
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

  model.type = 'gssm';                     % object type = generalized state space model
  model.tag  = 'GSSM_N1';                  % ID tag

  model.ffun      = @ffun;                   % file handle to FFUN
  model.hfun      = @hfun;                   % file handle to HFUN
  model.prior     = @prior;
  model.likelihood = @likelihood;            % file handle to LIKELIHOOD
  model.innovation = @innovation;            % file handle to INNOVATION
  model.linearize  = @linearize;              % file handle to LINEARIZE
  model.setparams  = @setparams;              % file handle to SETPARAMS

  model.statedim   = 1;                      %   state dimension
  model.obsdim     = 1;                      %   observation dimension
  model.paramdim   = 2;                      %   parameter dimension
  model.U1dim      = 1;                      %   exogenous control input 1 dimension
  model.U2dim      = 1;                      %   exogenous control input 2 dimension
  model.Vdim       = 1;                      %   process noise dimension
  model.Ndim       = 1;                      %   observation noise dimension

  Arg.type = 'gamma';
  Arg.dim = model.Vdim;
  Arg.alpha = 3;
  Arg.beta  = 0.5;
  model.pNoise = gennoiseds(Arg);   % process noise : Gamma(3,0.5) noise source

  Arg.type = 'gaussian';
  Arg.cov_type = 'full';
  Arg.dim = model.Ndim;
  Arg.mu = 0;
  Arg.cov  = 1e-5;
  Arg.cov_type = 'full';
  model.oNoise = gennoiseds(Arg);     % observation noise : zero mean white Gaussian noise, cov=0.2

  model.params = zeros(model.paramdim,1);

  model = setparams(model,[4e-2 0.5]);   % [omega phi]


%===============================================================================================
function model = setparams(model, params, index_vector)

  if (nargin==2)
    model.params = params(:);
  elseif (nargin==3)
    model.params(index_vector) = params(:);
  else
    error('[ setparams ] Incorrect number of input arguments.');
  end

%===============================================================================================
function new_state = ffun(model, state, V, U1)

  new_state      = 1 + sin(model.params(1)*pi.*U1) + model.params(2)*state;

  if ~isempty(V)
    new_state = new_state + V;
  end

%===============================================================================================
function observ = hfun(model, state, N, U2)

  [dim,nop] = size(state);

  observ = zeros(model.obsdim,nop);

  for k=1:nop,
    if (U2(k) <= 30),
       observ(k) = model.params(2)*state(:,k).^2;
    else
       observ = model.params(2)*state - 2;
    end
  end

  if ~isempty(N)
    observ = observ + N;
  end


%===============================================================================================
function tranprior = prior(model, nextstate, state, U1, pNoiseDS)

  X = nextstate - ffun(model, state, [], U1);

  tranprior = pNoiseDS.likelihood( pNoiseDS, X);


%===============================================================================================
function llh = likelihood(model, obs, state, U2, oNoiseDS)

  X = obs - hfun(model, state, [], U2);

  llh = oNoiseDS.likelihood( oNoiseDS, X);



%===============================================================================================
function innov = innovation(model, obs, observ)

  innov = obs - observ;


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
      out = model.params(2);

    case 'B'
      %%%========================================================
      %%%             Calculate B = df/dU1
      %%%========================================================
      out = [];

    case 'C'
      %%%========================================================
      %%%             Calculate C = dh/dx
      %%%========================================================
      if (U2 <= 30)
        out = 2*model.params(2)*state;
      else
        out = model.params(2);
      end

    case 'D'
      %%%========================================================
      %%%             Calculate D = dh/dU2
      %%%========================================================
      out = [];

    case 'G'
      %%%========================================================
      %%%             Calculate G = df/dv
      %%%========================================================
      out = 1;

    case 'H'
      %%%========================================================
      %%%             Calculate H = dh/dn
      %%%========================================================
      out = 1;

    case 'JFW'
      %%%========================================================
      %%%             Calculate  = dffun/dparameters
      %%%========================================================
      out = [cos(model.params(1)*pi*U1)*pi*U1 state];


    case 'JHW'
      %%%========================================================
      %%%             Calculate  = dhfun/dparameters
      %%%========================================================
      if (U2 <= 30)
        out = [0 state^2];
      else
        out = [0 state];
      end

    otherwise
      error('[ linearize ] Invalid model term requested!');

  end

  if (nargin==8), out = out(:,index_vector); end

  %--------------------------------------------------------------------------------------