% GSSM_SPEECH  Generalized state space model for single phoneme speech enhancement
%
% A single speech phoneme sampled at 8kHz is corrupted by additive colored (pink) noise.
% We use a simple linear autoregressive model (10th order) to model the generative model
% of the speech signal.
% We model the pink noise by a known 6th order linear autoregressive process driven by white Gaussian
% noise with known variance. The SNR of the noisy signal (y=clean+noise) is 0dB.
%
% The colored noise modeling (augmented state space model) is done according to the method proposed in:
% "Filtering of Colored Noise for Speech Enhancment and Coding", by J. D. Gibson, B. Koo and S. D. Gray,
% IEEE Transactions on Signal Processing, Vol. 39, No. 8, August 1991.
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


  load speech_data.mat noise_model noise_pnvar noisy clean;    % Loads colored noise model (LPC parameters) and process noise variance

  speech_taps = 10;
  speech_model = aryule(clean,speech_taps);
  speech_pnvar = var(filter(speech_model,1,clean));
  speech_model = -1*speech_model(2:end);
  noise_taps  = length(noise_model);     % number of noise filter taps

  %-- REQUIRED FIELDS

  model.type = 'gssm';                  % object type = generalized state space model
  model.tag  = 'GSSM_Speech_Colored_Noise_Linear';  % ID tag

  model.ffun       = @ffun;             % functionhandle to FFUN
  model.hfun       = @hfun;             % functionhandle to HFUN
  model.setparams  = @setparams;        % functionhandle to SETPARAMS

  model.statedim   = speech_taps + noise_taps;   % state dimension 10 for speech state + length of colored noise state
  model.obsdim     = 1;                 % observation dimension
  model.paramdim   = speech_taps + noise_taps;   % parameter dimension  (weights + colored noise parameters)
  model.U1dim      = 0;                 % exogenous control input 1 dimension
  model.U2dim      = 0;                 % exogenous control input 2 dimension
  model.Vdim       = 2;                 % process noise dimension  (augmented process noise needed for colored noise,
                                        % resulting in perfect measurment model with no explicit observation noise
  model.Ndim       = 0;                 % observation noise dimension (efective noise dimension is 0 for colored noise case)

  %-- SETUP NOISE DATA STRUCTURES

  Arg.type = 'gaussian';                % process noise source
  Arg.cov_type = 'full';
  Arg.dim = model.Vdim;
  Arg.mu = [0; 0];
  Arg.cov  = [speech_pnvar 0; 0 noise_pnvar];        % process noise variance
  model.pNoise = gennoiseds(Arg);       % generate process noise data structure : zero mean white Gaussian noise

  Arg.type = 'gaussian';
  Arg.dim = 0;
  Arg.mu = [];
  Arg.cov  = [];
  model.oNoise = gennoiseds(Arg);     % This observation noise model is actually only a dummy model, in that for the colored
                                      % noise case, the observation noise enters the state observation function implicitely.


  %-- OPTIONAL FIELDS

  model.noise_model = noise_model(:)';    % AR model for colored noise
  model.speech_model = speech_model(:)';
  model.speech_taps = speech_taps;
  model.noise_taps = noise_taps;

  %-- Call 'setparams' function once to make sure model parameters are correctly initialized

  model = setparams(model, [speech_model(:); noise_model(:)]);    % set/store the model parameters


%===============================================================================================
function model = setparams(model, params, index_vector)

  if (nargin==2)
    model.params = params(:);
  elseif (nargin==3)
    model.params(index_vector) = params(:);
  else
    error('[ setparams ] Incorrect number of input arguments.');
  end

  model.speech_model = model.params(1:model.speech_taps)';
  model.noise_model = model.params(model.speech_taps+1:end)';


%===============================================================================================
function new_state = ffun(model, state, V, U1)

  [dim,N] = size(state);

  speech_taps = model.speech_taps;

  new_state = zeros(dim,N);

  %-- SPEECH STATE UPDATE          -  linear AR
  new_state(1,:) = model.speech_model*state(1:speech_taps,:);
  new_state(2:speech_taps,:) = state(1:speech_taps-1,:);

  %-- COLORED NOISE STATE UPDATE   -  linear AR
  new_state(speech_taps+1,:) = model.noise_model*state(speech_taps+1:end,:);
  new_state(speech_taps+2:end,:) = state(speech_taps+1:end-1,:);

  if ~isempty(V)
    new_state(1,:) = new_state(1,:) + V(1,:);
    new_state(speech_taps+1,:) = new_state(speech_taps+1,:) + V(2,:);
  end

%===============================================================================================
function observ = hfun(model, state, N, U2)

  if isempty(N),

    observ = state(1,:) + state(model.speech_taps+1,:);

  else

    observ = state(1,:) + state(model.speech_taps+1,:) + N(1,:);

  end






