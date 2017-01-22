% GSSM  Template file for generalized state space model.
%
%   This template file is used to completely describe a system in a generalized
%   state space format useable by the ReBEL inference and estimation system.
%   This file must be copied, renamed and adapted to your specific problem. The
%   interface to each function should NOT BE CHANGED however.
%
%   The following main and subfunctions must be defined:
%
%   1) [VARARGOUT] = MODEL_INTERFACE(FUNC, VARARGIN) :  This function is the
%        main gateway function which is used to initialise the generalized state
%        space model data structure. This is done by calling the 'init'
%        subfunction. The user can extend this function to indirectly call other
%        subfunctions within this file if needed.
%
%   2) MODEL = INIT(INIT_ARGS) : This function generates and initializes a
%        generalized state space model (gssm) data structure which summarizes
%        all relevant information about the system. 'model' is a Matlab structure
%        that must contain the following fields (consistentency can be checked
%        with the 'consistent' function.
%       - type : String which contains the model type. Use 'gssm' for generalized
%                state space model.
%       - tag  : ID string which contains instance specific identification info.
%                Default value=''
%       - ffun : Function-handle to the state transition (state dynamics)
%                subfunction.
%       - hfun : Function-handle to the state observation subfunction.
%       - setparams : Function-handle to the setparams subfunction to update and
%                     possibly unpack the model parameters.
%
%       - prior : <<optional>> Function-handle to the state transition 'prior'
%                 function that calculates P(x(k)|x(k-1)). This must be defined
%                 if any of the particle filter family of estimators will be used
%                 on this model.
%       - likelihood : <<optional>> Function-handle to the observation likelihood
%                      function that calculates p(y(k)|x(k)) for a given
%                      realization of the state variable x and a particular
%                      observation instance y. This must be defined if any of the
%                      particle filter family of estimators will be used on this
%                      model.
%       - innovation : <<optional>> Function-handle to the innovation model
%                      function that calculates the difference between the output
%                      of the observation function (hfun) and the actual
%                      'real-world' measurement/observation of that signal. If
%                      this field is not defined, a generic innovation is used.
%       - linearize : <<optional>> Function-handle to the linearization
%                     subfunction. This is only needed if a linear Kalman filter
%                     (kf) or Extended Kalman Filter (ekf) will be used on this
%                     model. If no subfunction is defined, a default 'perturbation'
%                     based method of linearization will be used. This function
%                     does not need to be defined for the use of any of the
%                     Sigma-Point Kalman Filters (ukf, cdkf, srukf & srcdkf) or
%                     any of the Particle Filters (pf & sppf).
%
%       - statedim  : state dimension (this should be consistentent with the ffun
%                     and hfun subfunctions).
%       - obsdim    : observation dimension (this should be consistentent with
%                     the hfun subfunction).
%       - paramdim  : parameter dimension (number of free parameters in the
%                     system).
%       - U1dim     : dimension of exogenous input to ffun
%       - U2dim     : dimension of exogenous input to hfun
%       - pNoise    : process noise data structure  (this data structure is of
%                     type NoiseDS)
%       - oNoise    : observation noise data structure (this data structure is of
%                     type NoiseDS)
%       - params    : vector to hold all model parameters (must be of dimension
%                     [paramdim-by-1] )
%
%       - stateAngleCompIdxVec : <<optional>> Index vector idicating which (if
%                                any) of the state vector components are
%                                angular quantities (measured in radians) that
%                                has a discontinuety at +-Pi radians (this is
%                                needed by all SPKF based algorithms and derived
%                                hybrids)
%       - obsAngleCompIdxVec   : <<optional>> Index vector idicating which (if
%                                any) of the observation vector components are
%                                angular quantities (measured in radians) that
%                                has a discontinuety at +-Pi radians (this is
%                                needed by all SPKF based algorithms and
%                                derived hybrids)
%
%   3) MODEL = SETPARAMS(MODEL, PARAMS, IDXVECTOR) : This function unpacks a
%      column vector containing system parameters into specific forms needed by
%      FFUN, HFUN and possibly defined sub-functional objects. Both the
%      vectorized (packed) form of the parameters 'PARAMS' as well as the
%      unpacked forms are stored within the model data structure. 'IDXVECTOR' is
%      an optional argument which indicates which parameters should be updated.
%      This can be used to only modify a subset of the total system parameters.
%      'PARAMS' and 'IDXVECTOR' must have the same length.
%      Example :  model=setparams(model, [1 1 2 1 3]', [1 3 6:8])
%      << THIS SUBFUNCTION IS REQUIRED >>
%
%   4) NEW_STATE = FFUN(MODEL, STATE, V, U1) : State transition function which
%      takes as input the current state of the system 'STATE', a process noise
%      vector 'V', an exogenous control input 'U1', and a gssm data structure
%      'MODEL', and calculates the system state at the next discrete time instant,
%      'NEW_STATE'. This function implements the system dynamics.
%      << THIS SUBFUNCTION IS REQUIRED >>
%
%   5) OBSERV = HFUN(MODEL, STATE, N, U2) : State observation function which
%      takes as input the current state of the system 'STATE', an observation
%      noise vector 'N', an exogenous control input 'U2' and a gssm data
%      structure 'MODEL', and calculates the current observation vector of the
%      system, 'OBSERV'.
%      << THIS SUBFUNCTION IS REQUIRED >>
%
%   6) TRAN_PRIOR = PRIOR(MODEL, NEXT_STATE, STATE, U1, PNOISEDS) : Calculates
%      the transition prior p(next_state|state) = p(state(k)|state(k-1)) =
%      p(x(k)|x(k-1)) given a gssm data structure 'MODEL', realizations of the
%      system state at time k and k-1, 'NEXT_STATE' and 'STATE' and the
%      exogeneous inputs to the process model, U1. The process noise data
%      structure 'PNOISEDS' specifies which noise model should be used to
%      calculate the likelihood. If this is ommitted, the default model defined
%      process noise data structure 'model.pNoise' is used.
%      << THIS SUBFUNCTION IS OPTIONAL : Only required by particle filter
%         family of estimator >>
%
%   7) LLH = LIKELIHOOD(MODEL, OBS, STATE, U2, ONOISEDS) : Calculates the
%      likelihood of a 'real world' observation 'OBS' for a given realization
%      or instance of the state variable STATE. i.e. Calculates the value of
%      P(OBS|STATE). The measurement noise data structure 'ONOISEDS' specifies
%      which noise model should be used to calculate the likelihood. If this is
%      ommitted, the default model defined observation noise data structure
%      'model.pNoise' is used. 'U2' is the (optional) exogeneous input to the
%      state observation function 'hfun'.
%      << THIS SUBFUNCTION IS OPTIONAL : Only required by particle filter
%         family of estimator >>
%
%   8) INNOV = INNOVATION(MODEL, OBS, OBSERV) : Calculates the innovation signal
%      (difference) between the output of HFUN, i.e. OBSERV=HFUN(STATE) (the
%      predicted system observation) and an actual 'real world' observation OBS.
%      This function might be as simple as INNOV = OBS - OBSERV, which is the
%      default case, but can also be more complex for complex measurement
%      processes where for example multiple (possibly false) observations can be
%      observed for a given hidden ground truth.
%      << THIS SUBFUNCTION IS OPTIONAL : Only redefine if the default does not
%         reflect the true measurement process >>
%
%   9) OUT = LINEARIZE(MODEL, STATE, V, N, U1, U2, TERM, IDXVECTOR) generates a
%      linearized model of the nonlinear system described by the gssm data
%      structure MODEL at the current operating point, STATE, exogenous inputs U1
%      and U2. The linearized model is of the form:
%
%           state(k) = A*state(k-1) + B*u1(k-1) + G*v(k-1)
%               y(k) = C*state(k)   + D*u2(k)   + H*n(k)
%
%      for an arbitrary model defined by this GSSM file. The string TERM
%      specifies which of the model terms are returned, i.e.
%
%       A = linearize(model, state, v, n, u1, u2, 'A') or
%       H = linearize(model, state, v, n, u1, u2, 'H') etc.
%
%      TERM can be one of the following, 'A','B','C','D','G','H','JFW','JHW' ,
%      where 'JFW' and 'JHW' are the partial derivatives of FFUN and HFUN with
%      respect to the system parameters.
%
%      IDXVECTOR is an optional argument indicating which subset of the
%      independent vector should be used to calculate any specific derivative.
%      This will result in a Jacobian matrix with a reduced number of columns,
%      corresponding with the subvector as defined by 'IDXVECTOR'. The default
%      (when this argument is ommitted) is to use the full vector.
%
%      << THIS SUBFUNCTION IS OPTIONAL : Only required for Kalman and Extended
%      Kalman filters >>
%
%     See also
%     CONSIST, GENINFDS, GENNOISEDS
%
%   Copyright (c) Oregon Health & Science University (2006)
%
%   This file is part of the ReBEL Toolkit. The ReBEL Toolkit is available free for
%   academic use only (see included license file) and can be obtained from
%   http://choosh.csee.ogi.edu/rebel/.  Businesses wishing to obtain a copy of the
%   software should contact rebel@csee.ogi.edu for commercial licensing information.
%
%   See LICENSE.TXT (which should be part of the main toolkit distribution) for more
%   detail.


%===============================================================================================

function [varargout] = model_interface(func, varargin)

  switch func

    %--- Initialize GSSM data structure --------------------------------------------------------
    case 'init'
      model = init(varargin);
      error(consistent(model,'gssm'));              % check consistentency of initialized model
      varargout{1} = model;

    %--------------------------------------------------------------------------------------------
    otherwise

      error(['Function ''' func ''' not supported.']);

  end


%===============================================================================================
function model = init(init_args)

  model.type = 'gssm';                       % object type = generalized state space model

  model.tag  = '';                           % ID tag

  model.setparams = @setparams;              % function handle to SETPARAMS
  model.ffun      = @ffun;                   % function handle to FFUN
  model.hfun      = @hfun;                   % function handle to HFUN
  % model.prior = @prior;                    % function handle to PRIOR        (uncomment if 'prior' subfunction is defined)
  % model.likelihood = @likelihood;          % function handle to LIKELIHOOD   (uncomment if 'likelihood' subfunction is defined)
  % model.innovation = @innovation;          % function handle to INNOVATION   (uncomment if 'innovation' subfunction is defined)
  % model.linearize = @linearize;            % function handle to LINEARIZE    (uncomment if 'linearize' subfunction is defined)

  % model.stateAngleCompIdxVec :             % <<optional>> Index vector idicating which (if
                                             % any) of the state vector components are
                                             % angular quantities (measured in radians) that
                                             % has a discontinuety at +-Pi radians (this is
                                             % needed by all SPKF based algorithms and derived
                                             % hybrids)
  % model.obsAngleCompIdxVec                 % <<optional>> Same as model.stateAngleCompIdxVec but
                                             % used for the observation vector.

  %-- These dimensions have to be defined by the user

  model.statedim   = 0;                      % state dimension
  model.obsdim     = 0;                      % observation dimension
  model.paramdim   = 0;                      % total parameter dimension
  model.U1dim      = 0;                      % exogenous control input 1 dimension
  model.U2dim      = 0;                      % exogenous control input 2 dimension
  model.Vdim       = 0;                      % process noise dimension
  model.Ndim       = 0;                      % observation noise dimension


  %-- Setup process noise source

  Arg.type = 'gaussian';                     % noise source type : Gaussian
  Arg.cov_type = 'full';                      % Gaussian noise source cov_type (full covariance)
  Arg.tag = 'GSSM process noise source';     % Arbitrary ID tag (optional)
  Arg.dim = model.Vdim;                      % noise dimension
  Arg.mu = 0;                                % noise mean
  Arg.cov  = 1;                                % noise covariance
  model.pNoise = gennoiseds(Arg);            % generate noise source


  %-- Setup observation noise source

  Arg.type = 'gaussian';                     % noise source type : Gaussian
  Arg.cov_type = 'full';                      % Gaussian noise source cov_type (full covariance)
  Arg.tag = 'GSSM observation noise source'; % Arbitrary ID tag (optional)
  Arg.dim = model.Ndim;                      % noise dimension
  Arg.mu = 0;                                % noise mean
  Arg.cov  = 1;                                % noise covariance
  model.oNoise = gennoiseds(Arg);            % generate noise source


  model.params     = zeros(model.paramdim,1); %  setup parameter vector buffer

  model = setparams(model, zeros(model.paramdim,1));   % initialize model parameters and unpack if needed

  %-- The subsection of the init function below these comments should be used to define any other data structures, objects,
  %-- functions, etc.which is needed by the internal implementation of the FFUN, HFUN, LINEARIZE, ETC. functions. These data
  %-- structures should be saved within the GSSM 'model' data structure. The user can embed any other structures such as Netlab
  %-- neural networks, etc. in this section.


%===============================================================================================
function model = setparams(model, params, idxVector)

% Function to unpack a column vector containing system parameters into specific forms
% needed by FFUN, HFUN and possibly defined sub-functional objects. Both the vectorized (packed)
% form of the parameters as well as the unpacked forms are stored within the model data structure.
% INDEX_VECTOR is an optional argument which indicates which parameters should be updated. This can
% be used to only modify a subset of the total system parameters.
%
%   Example: model=setparams(model, [1 1 2 1 3]', [1 3 6:8]);

  switch nargin

  case 2
    %------------  Set all system parameters ---------------------------------------------------

    model.params = params;

    %-- Add unpack code here if needed -----
    %---------------------------------------


  case 3
    %------------  Set a subset of system parameters -------------------------------------------

    model.params(idxVector) = params;

    %-- Add unpack code if needed here.
    %---------------------------------------

  otherwise
    error('[ setparams ] Incorrect number of input arguments.');

  end


%===============================================================================================
function new_state = ffun(model, state, V, U1)

% FFUN  State transition function (system dynamics).
%
%   Generates the next state of the system NEW_STATE given
%   the current STATE, exogenous input U1 and process noise term V. If STATE, U1 and V are matrices
%   then FFUN is calculated for each column vector of these matrices, resulting in an equal number
%   of columns in NEW_STATE. MODEL is a GSSM derived data structure describing the system

%-- This function must be defined by the user!


%===============================================================================================
function  tranprior = prior(model, nextstate, state, U1, pNoiseDS)

% PRIOR  Transition prior function
%
%   Calculates P(nextstate|state). If you plan to run a particle filter on this mode, you should
%   define this.
%
%   INPUT
%         model          GSSM data structure
%         nextstate      state at time k
%         state          state at time k-1
%         U1             exogeneous input to FFUN at time k-1
%         pNoiseDS       (optional) process noise NoiseDS data structure to use for evaluation of
%                        transition prior. If this is ommitted, model.pNoise, is used.
%   OUTPUT
%         tranprior      p(x(k)|x(k-1))
%
%-- This function must be defined by the user!

%===============================================================================================
function observ = hfun(model, state, N, U2)

% HFUN  State observation function.
%
%   OBSERV = HFUN(MODEL, STATE, N, U2) generates the current possibly nonlinear observation of the
%   system state, OBSERV, given the current STATE, exogenous input U and observation noise term V.
%   If STATE, U2 and N are matrices then HFUN is calculated for each column vector of these matrices,
%   resulting in an equal number of columns in OBSERV. MODEL is a GSSM derived data structure describing
%   the system.

%-- This function must be defined by the user!


%===============================================================================================
function llh = likelihood(model, obs, state, U2, oNoiseDS)

% LIKELIHOOD  Observation likelihood function
%
% Function-handle to the observation likelihood function that calculates p(y|x) for a
% given realization of the state variable 'state' and a particular observation instance 'obs'.
%
%   i.e. Calculates the value of P(OBS|STATE) = P(y|x)
%
%   INPUT
%         model          GSSM data structure
%         obs            observation at time k
%         state          state at time k
%         U2             exogeneous input to HFUN at time k
%         oNoiseDS       (optional) measurement noise NoiseDS data structure to use for evaluation of
%                        transition prior. If this is ommitted, model.oNoise, is used.
%   OUTPUT
%         llh            p(y(k)|x(k))
%
%-- This function must be defined by the user!


%===============================================================================================
function INNOV = innovation(model, obs, observ)

% INNOVATION  Innovation model
%
%   INNOV = INNOVATION(MODEL, STATE, OBS, OBSERV) : Calculates the innovation signal (difference) between the
%   output of HFUN, i.e. OBSERV (the predicted system observation) and an actual 'real world' observation OBS.
%   This function might be as simple as INNOV = OBS - OBSERV, which is the default case, but can also be more
%   complex for complex measurement processes where for example multiple (possibly false) observations can be
%   observed for a given hidden ground truth.

%-- This function must be redefined by the user if the specific real world observation process dictates it


%===============================================================================================
function out = linearize(model, state, V, N, U1, U2, term, idxVector)

% LINEARIZE
%
%   OUT = LINEARIZE(MODEL, STATE, V, N, U1, U2, TERM, IDXVECTOR) returns a linearized model of the
%   form
%           state(k) = A*state(k-1) + B*u1(k-1) + G*v(k-1)
%               y(k) = C*state(k)   + D*u2(k)   + H*n(k)
%
%   for an arbitrary model defined by this GSSM file. The string TERM specifies which of the
%   model terms are returned, i.e.
%
%   A = linearize(model, state, v, n, u1, u2, 'A') or
%   O = linearize(model, state, v, n, u1, u2, 'H') etc.
%
%   TERM can be one of the following, 'A','B','C','D','G','H','JFW','JHW' , where 'JFW' and 'JHW'
%   are the partial derivatives of FFUN and HFUN with respect to the system parameters.
%
%   INDEX_VECTOR is an optional argument indicating which subset of the independent vector should be used to calculate
%   any specific derivative. This will result in a Jacobian matrix with a reduced number of columns,
%   corresponding with the subvector as defined by index_vector. The default (when this argument is ommitted)
%   is to use the full vector.
%
%   Generic perturbation based linearization subunits are provided. These can (and should) be replaced
%   by user defined analytical derivative code if available. If no linearization function is available
%   or is not needed, a call to this function should return an error message.
%

  nia = nargin;                      % number of input arguments
  if (nia < 7)
    error('[ linearize ] Not enough input arguments! ');
  end

  epsilon = 1e-8;                    % perturbation step size

  switch (term)

    case 'A'
      %%%========================================================
      %%%             Calculate A = dffun/dstate
      %%%========================================================
      if (nia==7), index_vector=[1:model.statedim]; end
      liv = length(index_vector);
      A  = zeros(model.statedim, liv);
      %%%---------- replace this section if needed --------------
      f1 = model.ffun(model,state,V,U1);
      for j=1:liv,
        s = state;
        k = index_vector(j);
        s(k) = s(k) + epsilon;
        f2 = model.ffun(model,s,V,U1);
        A(:,j) = (f2-f1)/epsilon;
      end
      %%%--------------------------------------------------------
      out = A;


    case 'B'
      %%%========================================================
      %%%             Calculate B = dffun/dU1
      %%%========================================================
      if (nia==7), index_vector=[1:model.U1dim]; end
      liv = length(index_vector);
      B = zeros(model.statedim, liv);
      %%%---------- replace this section if needed --------------
      f1 = model.ffun(model,state,V,U1);
      for j=1:liv,
        Utemp = U1;
        k = index_vector(j);
        Utemp(k) = Utemp(k) + epsilon;
        f2 = model.ffun(model,state,V,Utemp);
        B(:,j) = (f2-f1)/epsilon;
      end
      %%%--------------------------------------------------------
      out = B;


    case 'C'
      %%%========================================================
      %%%             Calculate C = dhfun/dx
      %%%========================================================
      if (nia==7), index_vector=[1:model.statedim]; end
      liv = length(index_vector);
      C = zeros(model.obsdim, liv);
      %%%---------- replace this section if needed --------------
      f3 = model.hfun(model,state,N,U2);
      for j=1:liv,
        s = state;
        k = index_vector(j);
        s(k) = s(k) + epsilon;
        f4 = model.hfun(model,s,N,U2);
        C(:,j) = (f4-f3)/epsilon;
      end
      %%%--------------------------------------------------------
      out = C;


    case 'D'
      %%%========================================================
      %%%             Calculate D = dhfun/dU2
      %%%========================================================
      if (nia==7), index_vector=[1:model.U2dim]; end
      liv = length(index_vector);
      D = zeros(model.obsdim, liv);
      %%%---------- replace this section if needed --------------
      f3 = model.hfun(model,state,N,U2);
      for j=1:liv,
        Utemp = U2;
        k = index_vector(j);
        Utemp(k) = Utemp(k) + epsilon;
        f4 = model.hfun(model,state,N,Utemp);
        D(:,j) = (f4-f3)/epsilon;
      end
      %%%--------------------------------------------------------
      out = D;


    case 'G'
      %%%========================================================
      %%%             Calculate G = dffun/dv
      %%%========================================================
      if (nia==7), index_vector=[1:model.Vdim]; end
      liv = length(index_vector);
      G = zeros(model.statedim, liv);
      %%%---------- replace this section if needed --------------
      f1 = model.ffun(model,state,V,U1);
      for j=1:liv,
        Vtemp = V;
        k = index_vector(j);
        Vtemp(k) = Vtemp(k) + epsilon;
        f5 = model.ffun(model,state,Vtemp,U1);
        G(:,j) = (f5-f1)/epsilon;
      end
      %%%--------------------------------------------------------
      out = G;


    case 'H'
      %%%========================================================
      %%%             Calculate H = dhfun/dn
      %%%========================================================
      if (nia==7), index_vector=[1:model.Ndim]; end
      liv = length(index_vector);
      H = zeros(model.obsdim, liv);
      %%%---------- replace this section if needed --------------
      f3 = model.hfun(model,state,N,U2);
      for j=1:liv,
        Ntemp = N;
        k = index_vector(j);
        Ntemp(k) = Ntemp(k) + epsilon;
        f6 = model.hfun(model,state,Ntemp,U2);
        H(:,j) = (f6-f3)/epsilon;
      end
      %%%--------------------------------------------------------
      out = H;


    case 'JFW'
      %%%========================================================
      %%%             Calculate  = dffun/dparameters
      %%%========================================================
      if (nia==7), index_vector=[1:model.paramdim]; end
      liv = length(index_vector);
      JFW = zeros(model.statedim, liv);
      %%%---------- replace this section if needed --------------
      f1 = model.ffun(model,state,V,U1);
      old_params = model.params;                         % save current model parameters
      for j=1:liv,
        params = old_params;
        k = index_vector(j);
        params(k) = params(k) + epsilon;
        model = setparams(model,params);
        f7 = model.ffun(model,state,V,U1);
        JFW(:,j) = (f7-f1)/epsilon;
      end
      %%%--------------------------------------------------------
      out = JFW;


    case 'JHW'
      %%%========================================================
      %%%             Calculate  = dhfun/dparameters
      %%%========================================================
      if (nia==7), index_vector=[1:model.paramdim]; end
      liv = length(index_vector);
      JHW = zeros(model.obsdim, liv);
      %%%---------- replace this section if needed --------------
      f3 = model.hfun(model,state,N,U2);
      old_params = model.params;                         % save current model parameters
      for j=1:liv,
        params = old_params;
        k = index_vector(j);
        params(k) = params(k) + epsilon;
        model = setparams(model,params);
        f8 = model.hfun(model,state,N,U2);
        JHW(:,j) = (f8-f3)/epsilon;
      end
      %%%--------------------------------------------------------
      out = JHW;

    otherwise
      error('[ linearize ] Invalid linearization term requested!');

  end

  %--------------------------------------------------------------------------------------
