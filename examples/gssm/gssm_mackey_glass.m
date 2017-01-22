% GSSM_MACKEY_GLASS  Generalized state space model for Mackey-Glass chaotic time series
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
%  tap length of 6 and a 6-4-1 MLP neural network with hyperbolic tangent activation functions
%  in the hidden layer and a linear output activation.
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

  load mg30_6-4-1_model.mat;            % load ReBEL neural network model from Matlab MAT file
                                        % This loads a NeuralNetDS data structure, i.e.
                                        % NeuralNet.type = 'NeuralNetDS'
                                        % NeuralNet.subtype = 'MLP'
                                        % NeuralNet.nodes   : MLP structure descriptor vector
                                        % NeuralNet.olType  : output layer type
                                        % NeuralNet.weights : neural network parameters
                                        % NeuralNet.pnVar   : inovation variance

  model.type = 'gssm';                  % object type = generalized state space model
  model.tag  = 'GSSM_Mackey-Glas-30';   % ID tag

  model.ffun      = @ffun;              % functionhandle to FFUN
  model.hfun      = @hfun;              % functionhandle to HFUN
  model.linearize = @linearize;         % functionhandle to LINEARIZE
  model.setparams = @setparams;         % functionhandle to SETPARAMS
  model.likelihood = @likelihood;       % functionhandle to LIKELIHOOD
  model.prior = @prior;                 % functionhandle to PRIOR function

  model.statedim   = 6;                 % state dimension
  model.obsdim     = 1;                 % observation dimension
  model.paramdim   = length(NeuralNet.weights);   % parameter dimension   [should equal 6*4+4 + 4*1 + 1 ]
  model.U1dim      = 0;                 % exogenous control input 1 dimension
  model.U2dim      = 0;                 % exogenous control input 2 dimension
  model.Vdim       = 1;                 % process noise dimension
  model.Ndim       = 1;                 % observation noise dimension


  Arg.type = 'gaussian';                % process noise source
  Arg.cov_type = 'full';
  Arg.dim = model.Vdim;
  Arg.mu = 0;
  Arg.cov  = NeuralNet.pnVar;                      % process noise variance
  model.pNoise = gennoiseds(Arg);       % process noise : zero mean white Gaussian noise

  Arg.type = 'gaussian';
  Arg.cov_type = 'full';
  Arg.dim = model.Ndim;
  Arg.mu = 0;
  Arg.cov  = 1;                           % This will be set in the main program, depending on the experiment
  model.oNoise = gennoiseds(Arg);

  model.params = zeros(model.paramdim,1);  % setup model parameter vector buffer (this is required by ReBEL)

  % Problem/model specific parameters are saved here to speed up subsequent access to these values. One can also
  % just save the whole neural network data structure, i.e. model.NeuralNetwork = NeuralNetwork, but this adds another
  % layer of dereferencing, which will slow down the code. This is up to the user to decide.

  model.nodes = NeuralNet.nodes;
  model.olType = NeuralNet.olType;
  model.trueWeights = NeuralNet.weights;

  % pre-allocate parameter buffers
  [model.W1, model.B1, model.W2, model.B2] = mlpunpack(model.nodes, model.trueWeights);

  % Generate NN parameter devectorizing indexes to allow for self-contained 'setparams' function. This
  % speeds up the code for parameter and joint estimation
  [model.idxW1, model.idxB1, model.idxW2, model.idxB2] = mlpindexgen(model.nodes);


  % Call setparam function (required)
  model = setparams(model, model.trueWeights, 1:model.paramdim);    % set/store the model parameters


%===============================================================================================
function model = setparams(model, params, paramIdxVec)

  switch nargin
   case 2
     model.params = params;
   case 3
     model.params(paramIdxVec) = params;
  end

  % Unpack ReBEL MLP Neural Net parameters 'inline'. This can also be accomplished with a call
  % to 'mlpunpack', but this way speeds up the code.

  tparams = model.params;

  model.W1(:) = tparams(model.idxW1);
  model.B1    = tparams(model.idxB1);
  model.W2(:) = tparams(model.idxW2);
  model.B2    = tparams(model.idxB2);


%===============================================================================================
function new_state = ffun(model, state, V, U1)

  nov = size(state,2);
  new_state = zeros(model.statedim,nov);

  % direct implementation of ReBEL MLP neural network call ... 'nnet2' can also be called, but this is faster
  new_state(1,:) = model.W2 * tanh(model.W1*state + cvecrep(model.B1,nov)) + cvecrep(model.B2,nov);
  new_state(2:end,:) = state(1:end-1,:);

  if ~isempty(V)
    new_state(1,:) = new_state(1,:) + V(1,:);
  end

%===============================================================================================
function observ = hfun(model, state, N, U2)

  observ = state(1,:);

  if ~isempty(N)
    observ = state(1,:) + N(1,:);
  end


%===============================================================================================
function tranprior = prior(model, nextstate, state, U1, pNoiseDS)

  X = nextstate - ffun(model, state, [], U1);

  tranprior = pNoiseDS.likelihood( pNoiseDS, X(1,:));


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

      A=zeros(model.statedim);                                          % quick init to zeros
      A(2:model.statedim,1:model.statedim-1) = eye(model.statedim-1);
      A(1,1:model.statedim) = mlpjacobian('dydx', model.olType, model.nodes, state, model.W1, model.B1, model.W2, model.B2);
      out = A;

    case 'B'
      %%%========================================================
      %%%             Calculate B = df/dU1
      %%%========================================================

      out = [];

    case 'C'
      %%%========================================================
      %%%             Calculate C = dh/dx
      %%%========================================================

      C = zeros(model.obsdim, model.statedim);
      C(1,1) = 1;
      out = C;

   case 'D'
      %%%========================================================
      %%%             Calculate D = dh/dU2
      %%%========================================================
      out = [];

    case 'G'
      %%%========================================================
      %%%             Calculate G = df/dv
      %%%========================================================
      G = zeros(model.statedim,1);
      G(1,1) = 1;
      out = G;

    case 'H'
      %%%========================================================
      %%%             Calculate H = dh/dn
      %%%========================================================
      H = zeros(model.obsdim,1);
      H(1,1) = 1;
      out = H;

    case 'JFW'
      %%%========================================================
      %%%             Calculate  = dffun/dparameters
      %%%========================================================
      JFW = zeros(model.statedim, model.paramdim);
      JFW(1,:) = mlpjacobian('dydw', model.olType, model.nodes, state, model.W1, model.B1, model.W2, model.B2);
      out = JFW;


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