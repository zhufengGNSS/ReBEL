function InferenceDS = geninfds(ArgDS)

% GENINFDS  Generate inference data structure from a generalized state space model and user defined inference parameters.
%
%   InferenceDS = geninfds(ArgDS)
%
%   This function generates a ReBEL inference data structure from a generalized state space model (gssm) data
%   structure. This inference data structure is used by all ReBEL estimation algorithms. All the parameters and
%   information needed to build the data structures are passed into the function by means of a user defined argument
%   data structure, ArgDS. The specifics of the input and output parameters are given
%   below.
%
%   INPUT ARGUMENT DATA STRUCTURE FIELDS  : ArgDS._______
%
%     .type               : (string)      Inference (estimation) type : 'state', 'parameter' or 'joint'
%     .tag                : (string)      Arbitrary user defined ID tag string
%     .model              : (gssm)        Generalized state space model descriptor (see gssm.m)
%
%     .paramParamIdxVec   : (r-vector)   <<OPTIONAL for parameter and joint estimation, not used by state estimation>>
%                                         Index vector containing indices of a subset of parameters which must be
%                                         estimated. Default (if ommitted) is to estimate all parameters.
%     .paramFunSelect     : (string)     <<OPTIONAL for parameter estimation, not used by state and joint estimation>>
%                                         String indicating which of the state-transition 'ffun' or state-observation
%                                         'hfun' functions should be used as the functional unit for parameter estimation.
%                                         The default value for this field (if ommitted) is 'both' which uses both functions,
%                                         i.e. obs=hfun(ffun(x)). One can also specify 'both-p' which uses both functions
%                                         in a parallel combination, i.e. obs = [ffun(x)]
%                                                                               [hfun(x)]
%     .paramFFunOutIdxVec : (r-vector)   <<OPTIONAL for parameter estimation, not used by state and joint estimation>>
%                                         Index vector containing indices of a subset of the output of FFUN which must be
%                                         used for parameter estimation observations. This argument is only used if
%                                         .paramFunSelect = 'both' or 'ffun'. The default is to use all FFUN outputs.
%     .paramHFunOutIdxVec : (r-vector)   <<OPTIONAL for parameter estimation, not used by state and joint estimation>>
%                                         Index vector containing indices of a subset of the output of HFUN which must be
%                                         used for parameter estimation observations. This argument is only used if
%                                         .paramFunSelect = 'both' or 'hfun'. The default is to use all HFUN outputs.
%
%
%   OUTPUT ARGUMENTS
%
%     InferenceDS         : inference data structure
%
%     See also
%     GSSM, GENNOISEDS, GENSYSNOISEDS, CONSIST
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

%========================================================================================================


%--- ERROR CHECKING -------------------------------------------------------------------------------------

if (nargin < 1)
    error(' [ geninfds ] Not enough inputs.');
end

if ~isstruct(ArgDS)
    error(' [ geninfds ] Input argument must be an argument data structure ArgDS.');
end

mf = checkstructfields(ArgDS,'type','model');  % Check existence of required data structure fields
if ~isempty(mf)
    error([' [ geninfds ] Argument data structure does not contain the following required fields : ' mf]);
end

if ~ischar(ArgDS.type)
    error(' [ geninfds ] ArgDS.type must be a string indicating the inference type, i.e. ''state'', ''parameter'' or ''joint''.');
end


%--- BUILD INFERENCE DATA STRUCTURE ----------------------------------------------------------------------

InferenceDS.type    = 'InferenceDS';   % data structure type
InferenceDS.inftype = ArgDS.type;      % inference type


if isfield(ArgDS,'tag'),
    InferenceDS.tag  = ArgDS.tag;                         % ID tag
else
    InferenceDS.tag  = '';
end

if ~isempty(consistent(ArgDS.model,'gssm'));                 % check for consistentency of GSSM data structure
    error([' [ geninfds ] There is an inconsistentency with the supplied GSSM data structure subfield. Please use ' ...
           'the CONSIST function on the GSSM model to determine the exact problem.']);
else
    model = ArgDS.model;
    InferenceDS.model = model;                            % embed GSSM model
end



switch ArgDS.type

    %---
    %--- STATE ESTIMATION  --------------------------------------------------------------------------------
    %---
    case 'state'

        %--- dimensions & other detail ---
        InferenceDS.statedim = model.statedim;                % state dimension
        InferenceDS.obsdim   = model.obsdim;                  % observation dimension
        InferenceDS.U1dim    = model.U1dim;                   % exogenous input 1 dimension
        InferenceDS.U2dim    = model.U2dim;                   % exogenous input 2 dimension
        InferenceDS.Vdim = model.pNoise.dim;                  % process noise dimension
        InferenceDS.Ndim = model.oNoise.dim;                  % observation noise dimension

        %--- function handles ---
        InferenceDS.ffun = @ffun_state;
        InferenceDS.hfun = @hfun_state;
        if isfield(model, 'prior'),
            InferenceDS.prior = @prior_state;
        end
        if isfield(model,'likelihood'),
            InferenceDS.likelihood = @likelihood_state;
        end
        if isfield(model,'innovation'),
            InferenceDS.innovation = @innovation_state;
        else
            InferenceDS.innovation = [];
        end
        if isfield(model,'linearize'),
            InferenceDS.linearize = @linearize_state;             % linearization function functionhandle
        else
            InferenceDS.linearize = @linearize_generic;           % generic (perturbation based linearization)
        end

        %--- other stuff ---
        % Index vectors indicating the presence of angular components in the state and observation vectors
        if isfield(model,'stateAngleCompIdxVec'),
            InferenceDS.stateAngleCompIdxVec = model.stateAngleCompIdxVec;
        else
            InferenceDS.stateAngleCompIdxVec = [];
        end
        if isfield(model,'obsAngleCompIdxVec'),
            InferenceDS.obsAngleCompIdxVec = model.obsAngleCompIdxVec;
        else
            InferenceDS.obsAngleCompIdxVec = [];
        end


    %---
    %--- PARAMETER ESTIMATION ----------------------------------------------------------------------------
    %---
    case 'parameter'

        % Check parameter index vector
        if ~isfield(ArgDS,'paramParamIdxVec')
            paramParamIdxVec = 1:model.paramdim;
        else
            paramParamIdxVec = ArgDS.paramParamIdxVec;
            % Check vector entries
            if ((max(paramParamIdxVec) > model.paramdim) | (min(paramParamIdxVec) < 1))
                error(' [ geninfds::parameter ] Parameter index vector has illegal entries');
            end
            % Check for duplicate index entries
            if checkdups(paramParamIdxVec)
                error(' [ geninfds::parameter ] Duplicate parameter index vector entries not allowed.');
            end
        end
        InferenceDS.paramParamIdxVec = paramParamIdxVec;                      % copy index vector in InferenceDS

        % Check parameter function select argument
        if ~isfield(ArgDS,'paramFunSelect')
            paramFunSelect = 'both';
        elseif stringmatch(ArgDS.paramFunSelect,{'ffun','hfun','both-p','both'})
            paramFunSelect = ArgDS.paramFunSelect;
        else
            error(' [ geninfds::parameter ] Unknown value used for paramFunSelect');
        end
        InferenceDS.paramFunSelect = paramFunSelect;                          % copy parameter function select string

        % Check ffun and hfun output index vectors
        if stringmatch(paramFunSelect,{'both-p','ffun'})
            if ~isfield(ArgDS,'paramFFunOutIdxVec')
                paramFFunOutIdxVec = 1:model.statedim;
            else
                paramFFunOutIdxVec = ArgDS.paramFFunOutIdxVec;
                % Check vector entries
                if ((max(paramFFunOutIdxVec) > model.statedim) | (min(paramFFunOutIdxVec) < 1))
                    error(' [ geninfds::parameter ] FFUN output index vector has illegal entries');
                end
                % Check for duplicate index entries
                if checkdups(paramFFunOutIdxVec)
                    error(' [ geninfds::parameter ] Duplicate FFUN output index vector entries not allowed.');
                end
            end
            InferenceDS.paramFFunOutIdxVec = paramFFunOutIdxVec;                  % copy ffun output
        end
        if stringmatch(paramFunSelect,{'both','hfun','both-p'})
            if ~isfield(ArgDS,'paramHFunOutIdxVec')
                paramHFunOutIdxVec = 1:model.obsdim;
            else
                paramHFunOutIdxVec = ArgDS.paramHFunOutIdxVec;
                % Check vector entries
                if ((max(paramHFunOutIdxVec) > model.obsdim) | (min(paramHFunOutIdxVec) < 1))
                    error(' [ geninfds::parameter ] HFUN output index vector has illegal entries');
                end
                % Check for duplicate index entries
                if checkdups(paramHFunIdxVec)
                    error(' [ geninfds::parameter ] Duplicate HFUN output index vector entries not allowed.');
                end
            end
            InferenceDS.paramHFunOutIdxVec = paramHFunOutIdxVec;                  % copy ffun output
        end



        %-- Setup rest of structure

        switch paramFunSelect

        %...................................................................................................................
        case 'both'

            %--- dimensions & other detail ---
            InferenceDS.statedim = length(paramParamIdxVec);                 % state dimension
            InferenceDS.obsdim = length(paramHFunOutIdxVec);                 % observation dimension
            InferenceDS.U1dim = 0;                                           % expgenous input 1 dimension
            InferenceDS.U2dim = model.U1dim + model.statedim + model.U2dim;  % exogenous input 2 dimension
            InferenceDS.Vdim = InferenceDS.statedim;                         % process noise dimension
            InferenceDS.Ndim = model.Vdim + model.Ndim;                      % observation noise dimension

            %--- functions ---
            InferenceDS.ffun = @ffun_parameter;                             % state transition function functionhandle
            InferenceDS.hfun = @hfun_parameter_both;                        % state observation function functionhandle
            if isfield(model,'linearize')
                InferenceDS.linearize = @linearize_parameter_both;          % linearization function functionhandle
            else
                InferenceDS.linearize = @linearize_generic;
            end
            InferenceDS.prior = @prior_parameter;
            if isfield(model,'likelihood'),
                InferenceDS.likelihood = @likelihood_parameter_both;
            end
            if isfield(model,'innovation'),
                InferenceDS.innovation = @innovation_state;
            else
                InferenceDS.innovation = [];
            end

            %--- copy/setup fixed linear model parameters ---
            InferenceDS.A = eye(InferenceDS.statedim);
            InferenceDS.B = [];
            InferenceDS.G = eye(InferenceDS.statedim);

            %--- other stuff ---
            % Index vectors indicating the presence of angular components in the state and observation vectors
            InferenceDS.stateAngleCompIdxVec = [];
            InferenceDS.obsAngleCompIdxVec = [];
            if isfield(model,'obsAngleCompIdxVec'),
               for k=1:length(model.obsAngleCompIdxVec),
                 idx = find(InferenceDS.paramHFunOutIdxVec == model.obsAngleCompIdxVec(k));
                 InferenceDS.obsAngleCompIdxVec = [InferenceDS.obsAngleCompIdxVec idx];
               end
            end


        %...................................................................................................................
        case 'both-p'

            %--- dimensions & other detail ---
            InferenceDS.statedim = length(paramParamIdxVec);                        % state dimension
            InferenceDS.obsdim = length(paramFFunOutIdxVec) + length(paramHFunOutIdxVec); % observation dimension
            InferenceDS.U1dim = 0;                                                  % exogenous input 1 dimension
            InferenceDS.U2dim = model.U1dim + 2*model.statedim + model.U2dim;       % exogenous input 2 dimension
            InferenceDS.Vdim = InferenceDS.statedim;                                % process noise dimension
            InferenceDS.Ndim = model.Vdim+model.Ndim;                               % observation noise dimension

            %--- functions ---
            InferenceDS.ffun = @ffun_parameter;                               % state transition function functionhandle
            InferenceDS.hfun = @hfun_parameter_bothp;                         % state observation function functionhandle
            if isfield(model,'linearize')
                InferenceDS.linearize = @linearize_parameter_bothp;           % linearization function functionhandle
            else
                InferenceDS.linearize = @linearize_generic;
            end
            InferenceDS.prior = @prior_parameter;
            if (isfield(model,'likelihood') & isfield(model,'prior')),
                InferenceDS.likelihood = @likelihood_parameter_bothp;
            end
            InferenceDS.innovation = [];

            %--- copy/setup fixed linear model parameters ---
            InferenceDS.A = eye(InferenceDS.statedim);
            InferenceDS.B = [];
            InferenceDS.G = eye(InferenceDS.statedim);

            %--- other stuff ---
            % Index vectors indicating the presence of angular components in the state and observation vectors
            InferenceDS.stateAngleCompIdxVec = [];
            InferenceDS.obsAngleCompIdxVec = [];
            if isfield(model,'stateAngleCompIdxVec'),
               for k=1:length(model.stateAngleCompIdxVec),
                 idx = find(InferenceDS.paramFFunOutIdxVec == model.stateAngleCompIdxVec(k));
                 InferenceDS.obsAngleCompIdxVec = [InferenceDS.obsAngleCompIdxVec idx];
               end
            end
            if isfield(model,'obsAngleCompIdxVec'),
               for k=1:length(model.obsAngleCompIdxVec),
                 idx = find(InferenceDS.paramHFunOutIdxVec == model.obsAngleCompIdxVec(k));
                 InferenceDS.obsAngleCompIdxVec = [InferenceDS.obsAngleCompIdxVec idx];
               end
            end


        %...................................................................................................................
        case 'ffun'

            %--- parameter dimensions ---
            InferenceDS.statedim = length(paramParamIdxVec);                  % state dimension
            InferenceDS.obsdim = length(paramFFunOutIdxVec);                  % observation dimension
            InferenceDS.U1dim = 0;                                            % exogenous input 1 dimension
            InferenceDS.U2dim = model.U1dim + model.statedim;                 % exogenous input 2 dimension
            InferenceDS.Vdim = InferenceDS.statedim;                          % process noise dimension
            InferenceDS.Ndim = model.Vdim;                                    % observation noise dimension

            %--- functions ---
            InferenceDS.ffun = @ffun_parameter;                               % state transition function functionhandle
            InferenceDS.hfun = @hfun_parameter_f;                             % state observation function functionhandle
            if isfield(model,'linearize')
                InferenceDS.linearize = @linearize_parameter_f;               % linearization function functionhandle
            else
                InferenceDS.linearize = @linearize_generic;
            end
            InferenceDS.prior = @prior_parameter;
            if isfield(model,'prior'),
                InferenceDS.likelihood = @likelihood_parameter_f;
            end
            InferenceDS.innovation = [];

            %--- copy/setup fixed linear model parameters ---
            InferenceDS.A         = eye(InferenceDS.statedim);
            InferenceDS.B         = [];
            InferenceDS.G         = eye(InferenceDS.statedim);

            %--- other stuff ---
            % Index vectors indicating the presence of angular components in the state and observation vectors
            InferenceDS.stateAngleCompIdxVec = [];
            InferenceDS.obsAngleCompIdxVec = [];
            if isfield(model,'stateAngleCompIdxVec'),
               for k=1:length(model.stateAngleCompIdxVec),
                 idx = find(InferenceDS.paramFFunOutIdxVec == model.stateAngleCompIdxVec(k));
                 InferenceDS.obsAngleCompIdxVec = [InferenceDS.obsAngleCompIdxVec idx];
               end
            end


        %...................................................................................................................
        case 'hfun'

            %--- parameter dimensions ---
            InferenceDS.statedim = length(paramParamIdxVec);                  % state dimension
            InferenceDS.obsdim = length(paramHFunOutIdxVec);                  % observation dimension
            InferenceDS.U1dim = 0;                                            % exogenous input 1 dimension
            InferenceDS.U2dim = model.U2dim + model.statedim;                 % exogenous input 2 dimension
            InferenceDS.Vdim = InferenceDS.statedim;                          % process noise dimension
            InferenceDS.Ndim = model.Ndim;                                    % observation noise dimension

            %--- functions ---
            InferenceDS.ffun      = @ffun_parameter;                          % state transition function functionhandle
            InferenceDS.hfun      = @hfun_parameter_h;                        % state observation function functionhandle
            if isfield(model,'linearize')
                InferenceDS.linearize = @linearize_parameter_h;               % linearization function functionhandle
            else
                InferenceDS.linearize = @linearize_generic;
            end
            InferenceDS.prior = @prior_parameter;
            if isfield(model,'likelihood'),
                InferenceDS.likelihood = @likelihood_parameter_h;
            end
            InferenceDS.innovation = [];

            %--- copy/setup fixed linear model parameters ---
            InferenceDS.A         = eye(InferenceDS.statedim);
            InferenceDS.B         = [];
            InferenceDS.G         = eye(InferenceDS.statedim);


            %--- other stuff ---
            % Index vectors indicating the presence of angular components in the state and observation vectors
            InferenceDS.stateAngleCompIdxVec = [];
            InferenceDS.obsAngleCompIdxVec = [];
            if isfield(model,'obsAngleCompIdxVec'),
               for k=1:length(model.obsAngleCompIdxVec),
                 idx = find(InferenceDS.paramHFunOutIdxVec == model.obsAngleCompIdxVec(k));
                 InferenceDS.obsAngleCompIdxVec = [InferenceDS.obsAngleCompIdxVec idx];
               end
            end


        otherwise
            error(' The only valid values for the funselect field are : ''both'' , ''ffun'' and ''hfun''.');

        end



    %---
    %--- JOINT ESTIMATION --------------------------------------------------------------------------------
    %---
    case 'joint'

        % Check parameter index vector
        if ~isfield(ArgDS,'paramParamIdxVec')
            paramParamIdxVec = 1:model.paramdim;
        else
            paramParamIdxVec = ArgDS.paramParamIdxVec;
            % Check vector entries
            if ((max(paramParamIdxVec) > model.paramdim) | (min(paramParamIdxVec) < 1))
                error(' [ geninfds::parameter ] Parameter index vector has illegal entries');
            end
            % Check for duplicate index entries
            if checkdups(paramParamIdxVec)
                error(' [ geninfds::parameter ] Duplicate parameter index vector entries not allowed.');
            end
        end
        InferenceDS.paramParamIdxVec = paramParamIdxVec;                  % copy index vector in InferenceDS

        %--- dimensions ---
        pdim = length(paramParamIdxVec);
        InferenceDS.paramParamIdxVec = paramParamIdxVec;                  % save index vector in InferenceDS
        InferenceDS.statedim  = model.statedim + pdim;                    % state dimension
        InferenceDS.obsdim    = model.obsdim;                             % observation dimension
        InferenceDS.U1dim     = model.U1dim;                              % exogenous input 1 dimension
        InferenceDS.U2dim     = model.U2dim;                              % exogenous input 2 dimension
        InferenceDS.Vdim      = model.Vdim + pdim;                        % process noise dimension
        InferenceDS.Ndim      = model.Ndim;                               % observation noise dimension

        %--- functions ---
        InferenceDS.ffun      = @ffun_joint;                              % state transition function functionhandle
        InferenceDS.hfun      = @hfun_joint;                              % state observation function functionhandle
        if isfield(model, 'prior'),
            InferenceDS.prior = @prior_joint;
        end
        if isfield(model, 'likelihood')
            InferenceDS.likelihood = @likelihood_joint;
        end
        if isfield(model, 'innovation')
            InferenceDS.innovation = @innovation_state;
        else
            InferenceDS.innovation = [];
        end
        if isfield(model, 'linearize')
            InferenceDS.linearize = @linearize_joint;                     % linearization function functionhandle
        else
            InferenceDS.linearize = @linearize_generic;
        end

        %--- other stuff ---
        % Index vectors indicating the presence of angular components in the state and observation vectors
        if isfield(model,'stateAngleCompIdxVec'),
            InferenceDS.stateAngleCompIdxVec = model.stateAngleCompIdxVec;
        else
            InferenceDS.stateAngleCompIdxVec = [];
        end
        if isfield(model,'obsAngleCompIdxVec'),
            InferenceDS.obsAngleCompIdxVec = model.obsAngleCompIdxVec;
        else
            InferenceDS.obsAngleCompIdxVec = [];
        end


%-----------------------------------------------------------------------------------------------------------------------------------------------
otherwise

  error([' [ geninfds ] Inference type ''' ArgDS.type ''' not supported.']);

end


% add misc. default data fields

InferenceDS.uType = 'TMU';             % Update type : Time-and-Measurement Update is the default
                                       % Other options are : TU - time update only or MU - measurement update only.



return;



%***********************************************************************************************
%***                                                                                         ***
%***                               SUB FUNCTION BLOCK                                        ***
%***                                                                                         ***
%***********************************************************************************************


%===============================================================================================
%================================== STATE ESTIMATION FUNCTIONS =================================

function new_state = ffun_state(InferenceDS, state, V, U1)

    %  FFUN_STATE  State transition function of meta system for state estimation
    %
    %    new_state = ffun_state(InferenceDS, state, V, U1)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %    OUTPUT
    %         new_state       : (c-vector) updated meta system state vector

    new_state = InferenceDS.model.ffun( InferenceDS.model, state, V, U1);

%-------------------------------------------------------------------------------------
function observ = hfun_state(InferenceDS, state, N, U2)

    %  HFUN_STATE  State observation function of meta system for state estimation
    %
    %    observ = hfun_state(InferenceDS, state, N, U2)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         N               : (c-vector) meta system observation noise vector
    %         U2              : (c-vector) meta system exogenous input 2
    %    OUTPUT
    %         observ          : (c-vector)  meta system observation vector

    observ = InferenceDS.model.hfun( InferenceDS.model, state, N, U2);

%-------------------------------------------------------------------------------------
function tran_prior = prior_state(InferenceDS, nextstate, state, U1, pNoiseDS)

    %  PRIOR_STATE  Calculates the transition prior probability P(x_k|x_(k-1))
    %
    %    tranprior = prior_state(InferenceDS, nextstate, state, pNoiseDS)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         nextstate       : (c-vector)  system state at time k
    %         state           : (c-vector)  system state at time k-1
    %         U1              : (c-vector) meta system exogenous input 1
    %         pNoiseDS        : (NoiseDS)   process noise data structure
    %    OUTPUT
    %         tranprior       : scalar probability P(x_k|x_(k-1))

    tran_prior = InferenceDS.model.prior( InferenceDS.model, nextstate, state, U1, pNoiseDS);


%-------------------------------------------------------------------------------------
function llh = likelihood_state(InferenceDS, obs, state, U2, oNoiseDS)

    %  LIKELIHOOD_STATE  Calculates the likelood of a real-world observation obs given
    %                    a realization of the predicted observation for a given state,
    %                    i.e. p(y|x) = p(obs|state)
    %
    %    llh = likelihood_state(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         state           : (c-vector)  meta system state vector
    %         U2              : (c-vector) meta system exogenous input 2
    %         oNoiseDS        : (NoiseDS)   observation noise data structure
    %    OUTPUT
    %         llh             : scalar  likelihood

    llh = InferenceDS.model.likelihood( InferenceDS.model, obs, state, U2, oNoiseDS);

%-------------------------------------------------------------------------------------
function innov = innovation_state(InferenceDS, obs, observ)

    %  INNOVATION_STATE  Calculates the innovation signal (difference) between the
    %   output of HFUN, i.e. OBSERV (the predicted system observation) and an actual
    %   'real world' observation OBS. This function might be as simple as
    %   INNOV = OBS - OBSERV, which is the default case, but can also be more
    %   complex for complex measurement processes where for example multiple (possibly false)
    %   observations can be observed for a given hidden ground truth.
    %
    %    innov = innovation_state(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         observ          : (c-vector)  meta system observation vector
    %    OUTPUT
    %         inov            : (c-vector) innovation sequence

    innov = InferenceDS.model.innovation( InferenceDS.model, obs, observ);

%-------------------------------------------------------------------------------------
function varargout = linearize_state(InferenceDS, state, V, N, U1, U2, varargin)

    %  LINEARIZE_STATE  Linearization function of meta system for state estimation
    %
    %    varargout = linearize_state(InferenceDS, state, V, N, U1, U2, varargin)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         N               : (c-vector) meta system observation noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %         U2              : (c-vector) meta system exogenous input 2
    %         varargin        : (strings) linearization terms wanted, e.g. 'A','B','G',....
    %    OUTPUT
    %         varargout       : (matrices) linearization terms corresponding with varargin strings

  nop = length(varargin);

  for k=1:nop,
    varargout{k} = InferenceDS.model.linearize( InferenceDS.model, state, V, N, U1, U2, varargin{k});
  end


%===========================================================================================================
%================================ PARAMETER ESTIMATION FUNCTIONS ===========================================

function new_state = ffun_parameter(InferenceDS, state, V, U1)

    %  FFUN_PARAMETER  State transition function of meta system for parameter estimation
    %
    %    new_state = ffun_parameter(InferenceDS, state, V, U1)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system system state vector
    %         V               : (c-vector) meta system process noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %    OUTPUT
    %         new_state       : (c-vector) updated meta system state vector
    %
    % Relationship between input arguments and external model (GSSM) variables
    %
    %   state -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof.
    %   U1    -> this is usually an empty matrix
    %   V     -> synthetic process noise (speeds up convergence)
    %

    new_state = state + V;

%-------------------------------------------------------------------------------------
function tran_prior = prior_parameter(InferenceDS, nextstate, state, U1, pNoiseDS)

    %  PRIOR_STATE  Calculates the transition prior probability P(x_k|x_(k-1))
    %
    %    tranprior = prior_parameter(InferenceDS, nextstate, state, pNoiseDS)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         nextstate       : (c-vector)  system state at time k
    %         state           : (c-vector)  system state at time k-1
    %         U1              : (c-vector) meta system exogenous input 1
    %         pNoiseDS        : (NoiseDS)   process noise data structure
    %    OUTPUT
    %         tranprior       : scalar probability P(x_k|x_(k-1))

    X = nextstate - state;

    tran_prior = pNoiseDS.likelihood( pNoiseDS, X);


%-------------------------------------------------------------------------------------
function observ = hfun_parameter_bothp(InferenceDS, state, N, U2)

    %  HFUN_PARAMETER_BOTHP  State observation function of meta system for parameter estimation using both ffun and hfun
    %                     from the underlying GSSM.
    %
    %    observ = hfun_parameter_bothp(InferenceDS, state, N, U2)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         N               : (c-vector) meta system observation noise vector
    %         U2              : (c-vector) meta system exogenous input 2
    %    OUTPUT
    %         observ          : (c-vector) meta system observation vector
    %
    % Relationship arguments and external model (GSSM) variables
    %
    %   state  -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof
    %   U2     -> [external_state(k-1) external_U1(k-1) external_state(k) external_U2(k)]'
    %   N      -> [external_process_noise(k-1) external_observation_noise(k)]'
    %   observ -> [external_state(k) external_observation(k)]'

    [dim,nov] = size(state);

    observ = zeros(InferenceDS.obsdim,nov);

    dimX  = InferenceDS.model.statedim;
%    dimO  = InferenceDS.model.obsdim;
    dimV  = InferenceDS.model.Vdim;
    dimN  = InferenceDS.model.Ndim;
    dimU1 = InferenceDS.model.U1dim;
%    dimU2 = InferenceDS.model.U2dim;

    ext_state_1     = U2(1:dimX,:);
    ext_proc_noise  = N(1:dimV,:);
    ext_U1          = U2(dimX+1:dimX+dimU1,:);
    ext_state_2     = U2(dimX+dimU1+1:dimX+dimU1+dimX,:);
    ext_obs_noise   = N(dimV+1:dimV+dimN,:);
    ext_U2          = U2(dimX+dimU1+dimX+1:end,:);

    ffun_idx = InferenceDS.paramFFunOutIdxVec;
    hfun_idx = InferenceDS.paramHFunOutIdxVec;

    dimF0 = length(ffun_idx);
    dimH0 = length(hfun_idx);

    % loop over all input vectors
    for k=1:nov,
        % set model parameter vector
        InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(:,k), InferenceDS.paramParamIdxVec);
        % FFUN part of observation
        FFunOut = InferenceDS.model.ffun( InferenceDS.model, ext_state_1(:,k), ext_proc_noise(:,k), ext_U1(:,k));
        % HFUN part of observation
        HFunOut = InferenceDS.model.hfun( InferenceDS.model, ext_state_2(:,k), ext_obs_noise(:,k), ext_U2(:,k));
        observ(1:dimF0,k) = FFunOut(ffun_idx);
        observ(dimF0+1:dimF0+dimH0,k) = HFunOut(hfun_idx);
    end

%-------------------------------------------------------------------------------------
function innov = innovation_parameter_bothp(InferenceDS, obs, observ)

    %  INNOVATION_PARAMETER_BOTHP  Calculates the innovation signal (difference) between the
    %   output of HFUN, i.e. OBSERV (the predicted system observation) and an actual
    %   'real world' observation OBS.
    %
    %    innov = innovation_parameter_bothp(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         observ          : (c-vector)  meta system observation vector
    %    OUTPUT
    %         inov            : (c-vector) innovation sequence

    [dim,nov] = size(observ);

    ffun_idx = InferenceDS.paramFFunOutIdxVec;

    dimF0 = length(ffun_idx);

    innov=zeros(InferenceDS.obsdim*nov);

    innov(1:dimF0,:) = obs(1:dimF0,:) - observ(1:dimF0,:);
    innov(dimF0+1:obsdim,:) = InferenceDS.model.innovation( InferenceDS.model, obs(dimF0+1:obsdim,:), ...
                                    observ(dimF0+1:obsdim,:));

%-------------------------------------------------------------------------------------
function llh = likelihood_parameter_bothp(InferenceDS, obs, state, U2, oNoiseDS)

    %  LIKELIHOOD_PARAMETER_BOTHP  Calculates the likelood of a real-world observation obs given
    %                           a realization of the predicted observation for a given state,
    %                           i.e. p(y|x) = p(obs|state)
    %
    %    llh = likelihood_parameter_bothp(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         state           : (c-vector)  meta system state vector
    %         U2              : (c-vector) meta system exogenous input 2
    %         oNoiseDS        : (NoiseDS)   observation noise data structure
    %    OUTPUT
    %         llh             : scalar  likelihood

    [dim,nov] = size(state);

    llh = zeros(1,nov);

    dimX  = InferenceDS.model.statedim;
%    dimO  = InferenceDS.model.obsdim;
    dimU1 = InferenceDS.model.U1dim;
%    dimU2 = InferenceDS.model.U2dim;

    ext_state_1     = U2(1:dimX,:);
%    ext_U1          = U2(dimX+1:dimX+dimU1,:);
    ext_state_2     = U2(dimX+dimU1+1:dimX+dimU1+dimX,:);
    ext_U2          = U2(dimX+dimU1+dimX+1:end,:);

    ffun_idx = InferenceDS.paramFFunOutIdxVec;
    hfun_idx = InferenceDS.paramHFunOutIdxVec;

    dimF0 = length(ffun_idx);
    dimH0 = length(hfun_idx);

    ext_nextstate = obs(1:dimF0,:);
    ext_obs = obs(dimF0+1:dimF0+dimH0,:);

    % loop over all input vectors
    for k=1:nov,

        % set model parameter vector
        InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(:,k), InferenceDS.paramParamIdxVec);

        % FFUN part of likelihood
        llh_f = InferenceDS.model.prior( InferenceDS.model, ext_nextstate(:,k), ext_state_1(:,k), oNoiseDS.noiseSources{1});

        % HFUN part of likelihood
        llh_h = InferenceDS.model.likelihood( InferenceDS.model, ext_obs(:,k), ext_state_2(:,k), ext_U2(:,k), oNoiseDS.noiseSources{2});

        llh(k) = llh_f * llh_h;       % we assume independence

    end


%-------------------------------------------------------------------------------------
function observ = hfun_parameter_f(InferenceDS, state, N, U2)

    %  HFUN_PARAMETER_F   State observation function of meta system for parameter estimation using only ffun
    %                     from the underlying GSSM.
    %
    %    observ = hfun_parameter_f(InferenceDS, state, N, U2)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         N               : (c-vector) meta system observation noise vector
    %         U2              : (c-vector) meta system exogenous input 2
    %    OUTPUT
    %         observ          : (c-vector) meta system observation vector
    %
    % Relationship between arguments and external model (GSSM) variables
    %
    %   state  -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof
    %   U2     -> [external_state(k-1) external_U1(k-1)]'
    %   N      -> [external_process_noise(k-1)]'
    %   observ -> [external_state(k)]'

    [dim,nov] = size(state);

    observ = zeros(InferenceDS.obsdim,nov);

    dimX  = InferenceDS.model.statedim;
    dimV  = InferenceDS.model.Vdim;
    dimU1 = InferenceDS.model.U1dim;

    ext_state_1     = U2(1:dimX,:);
    ext_proc_noise  = N(1:dimV,:);
    ext_U1          = U2(dimX+1:dimX+dimU1,:);

    ffun_idx = InferenceDS.paramFFunOutIdxVec;

%    dimF0 = length(ffun_idx);

    % loop over all input vectors
    for k=1:nov,
        % set model parameter vector
        InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(:,k), InferenceDS.paramParamIdxVec);
        FFunOut  = InferenceDS.model.ffun( InferenceDS.model, ext_state_1(:,k), ext_proc_noise(:,k), ext_U1(:,k));
        observ(:,k) = FFunOut(ffun_idx);
    end

%-------------------------------------------------------------------------------------
function observ = hfun_parameter_h(InferenceDS, state, N, U2)

    %  HFUN_PARAMETER_H   State observation function of meta system for parameter estimation using only hfun
    %                     from the underlying GSSM.
    %
    %    observ = hfun_parameter_h(InferenceDS, state, N, U2)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) system state vector
    %         N               : (c-vector) observation noise vector
    %         U2              : (c-vector) exogenous input 2
    %    OUTPUT
    %         observ          : (c-vector) observation vector
    %
    % Relationship between input arguments and external model (GSSM) variables
    %
    %   state  -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof
    %   U2     -> [external_state(k) external_U2(k)]'
    %   N      -> [external_observation_noise(k)]'
    %   observ -> [external_observation(k)]'

    [dim,nov] = size(state);

    observ = zeros(InferenceDS.obsdim,nov);

    dimX  = InferenceDS.model.statedim;
%    dimO  = InferenceDS.model.obsdim;
    dimN  = InferenceDS.model.Ndim;
    dimU2 = InferenceDS.model.U2dim;

    ext_state_2     = U2(1:dimX,:);
    ext_U2          = U2(dimX+1:dimX+dimU2,:);
    ext_obs_noise   = N(1:dimN,:);

    hfun_idx = InferenceDS.paramHFunOutIdxVec;

%    dimH0 = length(hfun_idx);

    % loop over all input vectors
    for k=1:nov,
       % set model parameter vector
       InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(:,k), InferenceDS.paramParamIdxVec);
       HFunOut = InferenceDS.model.hfun( InferenceDS.model, ext_state_2(:,k), ext_obs_noise(:,k), ext_U2(:,k));
       observ(:,k) = HFunOut(hfun_idx);
    end

%-------------------------------------------------------------------------------------
function observ = hfun_parameter_both(InferenceDS, state, N, U2)

    %  HFUN_PARAMETER_BOTH   State observation function of meta system for parameter estimation using the full system
    %                        dynamics of the underlying GSSM as observation, i.e. observ=hfun(ffun(x))
    %
    %
    %    observ = hfun_parameter_both(InferenceDS, state, N, U2)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) system state vector
    %         N               : (c-vector) observation noise vector
    %         U2              : (c-vector) exogenous input 2
    %    OUTPUT
    %         observ          : (c-vector) observation vector
    %
    % Relationship between input arguments and external model (GSSM) variables
    %
    %   state  -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof
    %   U2     -> [external_state(k-1) external_U1(k-1) external_U2(k)]'
    %   N      -> [external_observation_noise(k)]'
    %   observ -> [external_observation(k)]'

    [dim,nov] = size(state);

    observ = zeros(InferenceDS.obsdim,nov);

    dimX  = InferenceDS.model.statedim;
    dimV  = InferenceDS.model.Vdim;
    dimN  = InferenceDS.model.Ndim;
    dimU1 = InferenceDS.model.U1dim;
    dimU2 = InferenceDS.model.U2dim;

    ext_state_1     = U2(1:dimX,:);
    ext_U1          = U2(dimX+1:dimX+dimU1,:);
    ext_U2          = U2(dimX+dimU1+1:dimX+dimU1+dimU2,:);
    ext_obs_noise   = N(dimV+1:dimV+dimN,:);
    ext_proc_noise = N(1:dimV,:);
    
    hfun_idx = InferenceDS.paramHFunOutIdxVec;

%    dimH0 = length(hfun_idx);

    % loop over all input vectors
    for k=1:nov,
       % set model parameter vector
       InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(:,k), InferenceDS.paramParamIdxVec);
       % calculate X(k)=ffun(X(k-1))
       ext_state_2 = InferenceDS.model.ffun( InferenceDS.model, ext_state_1(:,k), ext_proc_noise(:,k), ext_U1(:,k));
       HFunOut = InferenceDS.model.hfun( InferenceDS.model, ext_state_2, ext_obs_noise(:,k), ext_U2(:,k));
       observ(:,k) = HFunOut(hfun_idx);
    end



%-------------------------------------------------------------------------------------
function llh = likelihood_parameter_f(InferenceDS, obs, state, U2, oNoiseDS)

    %  LIKELIHOOD_PARAMETER_F  Calculates the likelood of a real-world observation obs given
    %                           a realization of the predicted observation for a given state,
    %                           i.e. p(y|x) = p(obs|state)
    %
    %    llh = likelihood_parameter_f(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         state           : (c-vector)  meta system state vector
    %         U2              : (c-vector) meta system exogenous input 2
    %         oNoiseDS        : (NoiseDS)   observation noise data structure
    %    OUTPUT
    %         llh             : scalar  likelihood

    [dim,nov] = size(state);

    llh = zeros(1,nov);

    dimX  = InferenceDS.model.statedim;
    dimO  = InferenceDS.model.obsdim;
    dimU1 = InferenceDS.model.U1dim;

    ext_state_1     = U2(1:dimX,:);
    ext_U1          = U2(dimX+1:dimX+dimU1,:);

    ffun_idx = InferenceDS.paramFFunOutIdxVec;

    dimF0 = length(ffun_idx);

    ext_nextstate = obs(1:dimF0,:);

    % loop over all input vectors
    for k=1:nov,

        % set model parameter vector
        InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(:,k), InferenceDS.paramParamIdxVec);

        % FFUN part of likelihood
        llh(k) = InferenceDS.model.prior( InferenceDS.model, ext_nextstate(:,k), ext_state_1(:,k), ext_U1, oNoiseDS);

    end


%-------------------------------------------------------------------------------------
function llh = likelihood_parameter_h(InferenceDS, obs, state, U2, oNoiseDS)

    %  LIKELIHOOD_PARAMETER_H  Calculates the likelood of a real-world observation obs given
    %                           a realization of the predicted observation for a given state,
    %                           i.e. p(y|x) = p(obs|state)
    %
    %    llh = likelihood_parameter_h(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         state           : (c-vector)  meta system state vector
    %         U2              : (c-vector) meta system exogenous input 2
    %         oNoiseDS        : (NoiseDS)   observation noise data structure
    %    OUTPUT
    %         llh             : scalar  likelihood

    [dim,nov] = size(state);

    llh = zeros(1,nov);

    dimX  = InferenceDS.model.statedim;
%    dimO  = InferenceDS.model.obsdim;
%    dimU2 = InferenceDS.model.U2dim;

    ext_state_2     = U2(1:dimX,:);
    ext_U2          = U2(dimX+1:end,:);

%    hfun_idx = InferenceDS.paramHFunOutIdxVec;

%    dimH0 = length(hfun_idx);

    % loop over all input vectors
    for k=1:nov,

        % set model parameter vector
        InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(:,k), InferenceDS.paramParamIdxVec);

        llh(k) = InferenceDS.model.likelihood( InferenceDS.model, obs(:,k), ext_state_2(:,k), ext_U2(:,k), oNoiseDS);

    end


%-------------------------------------------------------------------------------------
function llh = likelihood_parameter_both(InferenceDS, obs, state, U2, oNoiseDS)

    %  LIKELIHOOD_PARAMETER_BOTH  Calculates the likelood of a real-world observation obs given
    %                             a realization of the predicted observation for a given state,
    %                             i.e. p(y|x) = p(obs|state)
    %
    %    llh = likelihood_parameter_both(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         state           : (c-vector)  meta system state vector
    %         U2              : (c-vector) meta system exogenous input 2
    %         oNoiseDS        : (NoiseDS)   observation noise data structure
    %    OUTPUT
    %         llh             : scalar  likelihood

    [dim,nov] = size(state);

    llh = zeros(1,nov);

    dimX  = InferenceDS.model.statedim;
%    dimO  = InferenceDS.model.obsdim;
    dimU1 = InferenceDS.model.U1dim;
%    dimU2 = InferenceDS.model.U2dim;

    ext_state_1     = U2(1:dimX,:);
    ext_U1          = U2(dimX+1:dimX+dimU1,:);
    ext_U2          = U2(dimX+dimU1+1:end,:);

%    hfun_idx = InferenceDS.paramHFunOutIdxVec;

%    dimH0 = length(hfun_idx);

        % loop over all input vectors
    for k=1:nov,

        % set model parameter vector
        InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(:,k), InferenceDS.paramParamIdxVec);

        ext_state_2 = InferenceDS.model.ffun( InferenceDS.model, ext_state_1(:,k), [], ext_U1(:,k));

        llh(k) = InferenceDS.model.likelihood( InferenceDS.model, obs(:,k), ext_state_2, ext_U2(:,k), oNoiseDS);

    end

%--------------------------------------------------------------------------------------
function varargout = linearize_parameter_both(InferenceDS, state, V, N, U1, U2, varargin)

    %  LINEARIZE_PARAMETER_BOTH  Linearization function of meta system for parameter estimation using both ffun
    %                            and hfun from the underlying GSSM in a
    %                            cascading (i.e. y=hfun(ffun(state,U1,V),U2,N))
    %
    %    varargout = linearize_parameter_both(InferenceDS, state, V, N, U1, U2, varargin)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         N               : (c-vector) meta system observation noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %         U2              : (c-vector) meta system exogenous input 2
    %         varargin        : (strings) linearization terms wanted, e.g. 'A','B','G',....
    %    OUTPUT
    %         varargout       : (matrices) linearization terms corresponding with varargin strings
    %
    % Relationship between input arguments and external model (GSSM) variables
    %
    %   state -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof
    %   U1    -> this is usually an empty matrix
    %   U2     -> [external_state(k-1) external_U1(k-1) external_U2(k)]'
    %   V     -> synthetic process noise (speeds up convergence)
    %   N     -> [external_process_noise(k-1) external_observation_noise(k)]'


    % Setup temporary model to use for linearization purposes
    model = InferenceDS.model;                                                      % copy existing model
    if ~isempty(state),
        model = model.setparams( model, state, InferenceDS.paramParamIdxVec);   % set parameters acording to state variable
    end

    dimX  = model.statedim;
%    dimO  = model.obsdim;
    dimV  = model.Vdim;
    dimN  = model.Ndim;
    dimU1 = model.U1dim;
%    dimU2 = model.U2dim;

    ext_state_1     = U2(1:dimX);
    ext_proc_noise  = N(1:dimV);
    ext_U1          = U2(dimX+1:dimX+dimU1);
    ext_obs_noise   = N(dimV+1:dimV+dimN);
    ext_U2          = U2(dimX+dimU1+1:end);

%    ffun_idx = InferenceDS.paramFFunOutIdxVec;
    hfun_idx = InferenceDS.paramHFunOutIdxVec;

%    dimF0 = length(ffun_idx);
    dimH0 = length(hfun_idx);

    for k=1:length(varargin)

        switch varargin{k}

        %--- A = dffun/dstate
        case 'A'
            varargout{k} = InferenceDS.A;

        %--- B = dffun/dU1
        case 'B'
            varargout{k} = InferenceDS.B;

        %--- G = dffun/dv
        case 'G'
            varargout{k} = InferenceDS.G;

        %--- C = dhfun/dstate
        case 'C'
            C = zeros(InferenceDS.obsdim, InferenceDS.statedim);
            ext_state_2 = model.ffun( model, ext_state_1, ext_proc_noise, ext_U1);
            extC = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'C');
            extJFW = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'JFW', InferenceDS.paramParamIdxVec);
            extJHW = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'JHW', InferenceDS.paramParamIdxVec);
            Ctemp = extC*extJFW + extJHW;
            C(1:dimH0,:) = Ctemp(hfun_idx,:);
            varargout{k} = C;

        %--- D = dhfun/dU2
        case 'D'
            D = zeros(InferenceDS.obsdim, InferenceDS.U2dim);
            ext_state_2 = model.ffun( model, ext_state_1, ext_proc_noise, ext_U1);
            extA = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'A');
            extB = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'B');
            extC = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'C');
            extD = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'D');
            tempCA = extC*extA;
            tempCB = extC*extB;
            D(1:dimH0,1:dimX) = tempCA(hfun_idx,:);
            D(1:dimH0,dimX+1:dimX+dimU1) = tempCB(hfun_idx,:);
            D(1:dimH0,dimX+dimU1+1:end) = extD(hfun_idx,:);
            varargout{k} = D;

        %--- H = dhfun/dn
        case 'H'
            H = zeros(InferenceDS.obsdim, InferenceDS.Ndim);
            ext_state_2 = model.ffun( model, ext_state_1, ext_proc_noise, ext_U1);
            extC = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'C');
            extG = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'G');
            extH = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'H');
            tempCG = extC*extG;
            H(1:dimH0,1:dimV) = tempCG(hfun_idx,:);
            H(1:dimH0,dimV+1:end) = extH(hfun_idx,:);
            varargout{k} = H;

        %----
        otherwise
            error('[ InferenceDS.linearize ] Unknown linearization term.');

        end

    end


%--------------------------------------------------------------------------------------
function varargout = linearize_parameter_bothp(InferenceDS, state, V, N, U1, U2, varargin)

    %  LINEARIZE_PARAMETER_BOTHP  Linearization function of meta system for parameter estimation using both ffun
    %                          and hfun from the underlying GSSM.
    %
    %    varargout = linearize_parameter_bothp(InferenceDS, state, V, N, U1, U2, varargin)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         N               : (c-vector) meta system observation noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %         U2              : (c-vector) meta system exogenous input 2
    %         varargin        : (strings) linearization terms wanted, e.g. 'A','B','G',....
    %    OUTPUT
    %         varargout       : (matrices) linearization terms corresponding with varargin strings
    %
    % Relationship between input arguments and external model (GSSM) variables
    %
    %   state -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof
    %   U1    -> this is usually an empty matrix
    %   U2    -> [external_state(k-1) external_U1(k-1) external_state(k) external_U2(k)]'
    %   V     -> synthetic process noise (speeds up convergence)
    %   N     -> [external_process_noise(k-1) external_observation_noise(k)]'


    % Setup temporary model to use for linearization purposes
    model = InferenceDS.model;                                                      % copy existing model
    if ~isempty(state),
        model = model.setparams( model, state, InferenceDS.paramParamIdxVec);   % set parameters according to state variable
    end

    dimX  = model.statedim;
%    dimO  = model.obsdim;
    dimV  = model.Vdim;
    dimN  = model.Ndim;
    dimU1 = model.U1dim;
%    dimU2 = model.U2dim;

    ext_state_1     = U2(1:dimX);
    ext_proc_noise  = N(1:dimV);
    ext_U1          = U2(dimX+1:dimX+dimU1);
    ext_state_2     = U2(dimX+dimU1+1:dimX+dimU1+dimX);
    ext_obs_noise   = N(dimV+1:dimV+dimN);
    ext_U2          = U2(dimX+dimU1+dimX+1:end);

    ffun_idx = InferenceDS.paramFFunOutIdxVec;
    hfun_idx = InferenceDS.paramHFunOutIdxVec;

    dimF0 = length(ffun_idx);
    dimH0 = length(hfun_idx);

    for k=1:length(varargin)

        switch varargin{k}

        %--- A = dffun/dstate
        case 'A'
            varargout{k} = InferenceDS.A;

        %--- B = dffun/dU1
        case 'B'
            varargout{k} = InferenceDS.B;

        %--- G = dffun/dv
        case 'G'
            varargout{k} = InferenceDS.G;

        %--- C = dhfun/dstate
        case 'C'
            C = zeros(InferenceDS.obsdim, InferenceDS.statedim);
            extJFW = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'JFW', InferenceDS.paramParamIdxVec);
            extJHW = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'JHW', InferenceDS.paramParamIdxVec);
            C(1:dimF0,:) = extJFW(ffun_idx,:);
            C(dimF0+1:dimF0+dimH0,:) = extJHW(hfun_idx,:);
            varargout{k} = C;

        %--- D = dhfun/dU2
        case 'D'
            D = zeros(InferenceDS.obsdim, InferenceDS.U2dim);
            extA = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'A');
            extB = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'B');
            extC = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'C');
            extD = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'D');
            D(1:dimF0,1:dimX) = extA(ffun_idx,:);
            D(1:dimF0,dimX+1:dimX+dimU1) = extB(ffun_idx,:);
            D(dimF0+1:dimF0+dimH0,dimX+dimU1+1:dimX+dimU1+dimX) = extC(hfun_idx,:);
            D(dimF0+1:dimF0+dimH0,dimX+dimU1+dimX+1:end) = extD(hfun_idx,:);
            varargout{k} = D;

        %--- H = dhfun/dn
        case 'H'
            H = zeros(InferenceDS.obsdim, InferenceDS.Ndim);
            extG = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'G');
            extH = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'H');
            H(1:dimF0,1:dimV) = extG(ffun_idx,:);
            H(dimF0+1:dimF0+dimH0,dimV+1:end) = extH(hfun_idx,:);
            varargout{k} = H;

        %----
        otherwise
            error('[ InferenceDS.linearize ] Unknown linearization term.');

        end

    end

%-------------------------------------------------------------------------------------
function varargout = linearize_parameter_f(InferenceDS, state, V, N, U1, U2, varargin)

    %  LINEARIZE_PARAMETER_F  Linearization function of meta system for parameter estimation using only
    %                         ffun from the underlying GSSM.
    %
    %    varargout = linearize_parameter_f(InferenceDS, state, V, N, U1, U2, varargin)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         N               : (c-vector) meta system observation noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %         U2              : (c-vector) meta system exogenous input 2
    %         varargin        : (strings) linearization terms wanted, e.g. 'A','B','G',....
    %    OUTPUT
    %         varargout       : (matrices) linearization terms corresponding with varargin strings
    %
    % Relationship between input arguments and external model (GSSM) variables
    %
    %   state -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof
    %   U1    -> this is usually an empty matrix
    %   U2    -> [external_state(k-1) external_U1(k-1)]'
    %   V     -> synthetic process noise (speeds up convergence)
    %   N     -> [external_process_noise(k-1)]'


    % Setup temporary model to use for linearization purposes
    model = InferenceDS.model;                                                        % copy existing model
    if ~isempty(state),
        model = model.setparams( model, state, InferenceDS.paramParamIdxVec);   % set parameters acording to state variable
    end

    dimX  = model.statedim;
    dimV  = model.Vdim;
    dimU1 = model.U1dim;

    if isempty(U2)
        ext_state_1     = [];
        ext_U1          = [];
    else
        ext_state_1     = U2(1:dimX);
        ext_U1          = U2(dimX+1:dimX+dimU1);
    end
    if isempty(N),
        ext_proc_noise  = [];
    else
        ext_proc_noise  = N(1:dimV);
    end

    ffun_idx = InferenceDS.paramFFunOutIdxVec;

    dimF0 = length(ffun_idx);


    for k=1:length(varargin)

        switch varargin{k}

        %--- A = dffun/dstate
        case 'A'
            varargout{k} = InferenceDS.A;

        %--- B = dffun/dU1
        case 'B'
            varargout{k} = InferenceDS.B;

        %--- G = dffun/dv
        case 'G'
            varargout{k} = InferenceDS.G;

        %--- C = dhfun/dstate
        case 'C'
            C = zeros(InferenceDS.obsdim, InferenceDS.statedim);
            extJFW = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'JFW', InferenceDS.paramParamIdxVec);
            C = extJFW(ffun_idx,:);
            varargout{k} = C;

        %--- D = dhfun/dU2
        case 'D'
            D = zeros(InferenceDS.obsdim, InferenceDS.U2dim);
            extA = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'A');
            extB = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'B');
            D(1:dimF0,1:dimX) = extA(ffun_idx,:);
            D(1:dimF0,dimX+1:dimX+dimU1) = extB(ffun_idx,:);
            varargout{k} = D;

        %--- H = dhfun/dn
        case 'H'
            H = zeros(InferenceDS.obsdim, InferenceDS.Ndim);
            extG = model.linearize( model, ext_state_1, ext_proc_noise, [], ext_U1, [], 'G');
            H(1:dimF0,1:dimV) = extG(ffun_idx,:);
            varargout{k} = H;

        %---
        otherwise
            error('[ InferenceDS.linearize ] Unknown linearization term.');

        end

    end


%-------------------------------------------------------------------------------------
function varargout = linearize_parameter_h(InferenceDS, state, V, N, U1, U2, varargin)

    %  LINEARIZE_PARAMETER_H  Linearization function of meta system for parameter estimation using
    %                         only hfun from the underlying GSSM.
    %
    %    varargout = linearize_parameter_h(InferenceDS, state, V, N, U1, U2, varargin)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         N               : (c-vector) meta system observation noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %         U2              : (c-vector) meta system exogenous input 2
    %         varargin        : (strings) linearization terms wanted, e.g. 'A','B','G',....
    %    OUTPUT
    %         varargout       : (matrices) linearization terms corresponding with varargin strings
    %
    % Relationship between input arguments and external model (GSSM) variables
    %
    %   state -> external model parameters or a subset (specified by InferenceDS.paramParamIdxVec) thereof
    %   U1    -> this is usually an empty matrix
    %   U2    -> [external_state(k) external_U2(k)]'
    %   V     -> synthetic process noise (speeds up convergence)
    %   N     -> [external_observation_noise(k)]'


    % Setup temporary model to use for linearization purposes
    model = InferenceDS.model;                                                     % copy existing model
    model = model.setparams( model, state, InferenceDS.paramParamIdxVec);       % set parameters acording to state variable

    dimX  = model.statedim;
%    dimO  = model.obsdim;
    dimN  = model.Ndim;
    dimU2 = model.U2dim;

    ext_state_2     = U2(1:dimX);
    ext_obs_noise   = N(1:dimN);
    ext_U2          = U2(dimX+1:dimX+dimU2);

    hfun_idx = InferenceDS.paramHFunOutIdxVec;

    dimH0 = length(hfun_idx);


    for k=1:length(varargin),

        switch varargin{k}

        %--- A = dffun/dstate
        case 'A'
            varargout{k} = InferenceDS.A;

        %--- B = dffun/dU1
        case 'B'
            varargout{k} = InferenceDS.B;

        %--- G = dffun/dv
        case 'G'
            varargout{k} = InferenceDS.G;


        %--- C = dhfun/dstate
        case 'C'
%            C = zeros(InferenceDS.obsdim, InferenceDS.statedim);
            extJHW = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'JHW', InferenceDS.paramParamIdxVec);
            varargout{k} = extJHW(hfun_idx,:);


        %--- D = dhfun/dU2
        case 'D'
            D = zeros(InferenceDS.obsdim, InferenceDS.U2dim);
            extC = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'C');
            extD = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'D');
            D(1:dimH0,1:dimX) = extC(hfun_idx,:);
            D(1:dimH0,dimX+1:dimX+dimU2) = extD(hfun_idx,:);
            varargout{k} = D;


        %--- H = dhfun/dn
        case 'H'
            H = zeros(InferenceDS.obsdim, InferenceDS.Ndim);
            extH = model.linearize( model, ext_state_2, [], ext_obs_noise, [], ext_U2, 'H');
            H(1:dimH0,1:dimN) = extH(hfun_idx,:);
            varargout{k} = H;

        %---
        otherwise
            error('[ InferenceDS.linearize ] Unknown linearization term.');

        end

    end


%===================================================================================================
%================================= JOINT ESTIMATION FUNCTIONS ======================================

function new_state = ffun_joint(InferenceDS, state, V, U1)

    %  FFUN_JOINT  State transition function of meta system for joint estimation
    %
    %    new_state = ffun_joint(InferenceDS, state, V, U1)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %    OUTPUT
    %         new_state       : (c-vector) updated meta system state vector

    [dim,nov] = size(state);

    new_state = zeros(dim,nov);

    dimX  = InferenceDS.model.statedim;
    dimV  = InferenceDS.model.Vdim;

    for k=1:nov,
      InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(dimX+1:end,k), InferenceDS.paramParamIdxVec);    % set model parameter vector
      new_state(1:dimX,k) = InferenceDS.model.ffun( InferenceDS.model, state(1:dimX,k), V(1:dimV,k), U1(:,k));
      new_state(dimX+1:end,k) = state(dimX+1:end,k) + V(dimV+1:end,k);
    end

%-------------------------------------------------------------------------------------
function tran_prior = prior_joint(InferenceDS, nextstate, state, U1, pNoiseDS)

    %  PRIOR_JOINT  Calculates the transition prior probability P(x_k|x_(k-1))
    %
    %    tranprior = prior_joint(InferenceDS, nextstate, state, pNoiseDS)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         nextstate       : (c-vector)  system state at time k
    %         state           : (c-vector)  system state at time k-1
    %         U1              : (c-vector) meta system exogenous input 1
    %         pNoiseDS        : (NoiseDS)   process noise data structure
    %    OUTPUT
    %         tranprior       : scalar probability P(x_k|x_(k-1))


    [dim,nov] = size(state);

    tran_prior = zeros(1,nov);

    dimX  = InferenceDS.model.statedim;

    paramState = state(dimX+1:end,:);
    paramNextState = nextstate(dimX+1:end,:);

    stateState = state(1:dimX,:);
    stateNextState = nextstate(1:dimX,:);

    dX = paramNextState - paramState;

    tran_prior = pNoiseDS.likelihood( pNoiseDS, dX, pNoiseDS.N);

    for k=1:nov,
        InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, paramState(:,k), InferenceDS.paramParamIdxVec); % set model parameter vector
        tran_prior(k) = tran_prior(k) * InferenceDS.model.prior( InferenceDS.model, stateNextState(:,k), stateState(:,k), U1(:,k), pNoiseDS.noiseSources{1});
    end

%-------------------------------------------------------------------------------------
function observ = hfun_joint(InferenceDS, state, N, U2)

    %  HFUN_JOINT  State observation function of meta system for joint estimation
    %
    %    observ = hfun_joint(InferenceDS, state, N, U2)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         N               : (c-vector) meta system observation noise vector
    %         U2              : (c-vector) meta system exogenous input 2
    %    OUTPUT
    %         observ          : (c-vector)  meta system observation vector

    [dim,nov] = size(state);

    observ = zeros(InferenceDS.obsdim,nov);

    dimX  = InferenceDS.model.statedim;
%    dimV  = InferenceDS.model.Vdim;

    for k=1:nov,
      InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(dimX+1:end,k), InferenceDS.paramParamIdxVec);    % set model parameter vector
      observ(:,k) = InferenceDS.model.hfun( InferenceDS.model, state(1:dimX,k), N(:,k), U2(:,k));
    end

%-------------------------------------------------------------------------------------
function llh = likelihood_joint(InferenceDS, obs, state, U2, oNoiseDS)

    %  LIKELIHOOD_JOINT  Calculates the likelood of a real-world observation obs given
    %                           a realization of the predicted observation for a given state,
    %                           i.e. p(y|x) = p(obs|state)
    %
    %    llh = likelihood_joint(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         state           : (c-vector)  meta system state vector
    %         U2              : (c-vector) meta system exogenous input 2
    %         oNoiseDS        : (NoiseDS)   observation noise data structure
    %    OUTPUT
    %         llh             : scalar  likelihood

    [dim,nov] = size(state);

    llh = zeros(1,nov);

    dimX  = InferenceDS.model.statedim;

    for k=1:nov,
      InferenceDS.model = InferenceDS.model.setparams( InferenceDS.model, state(dimX+1:end,k), InferenceDS.paramParamIdxVec); ...
      % set model parameter vector
      llh(k) = InferenceDS.model.likelihood( InferenceDS.model, obs(:,k), state(1:dimX,k), U2(:,k), oNoiseDS);
    end


%-------------------------------------------------------------------------------------
function varargout = linearize_joint(InferenceDS, state, V, N, U1, U2, varargin)

    %  LINEARIZE_JOINT  Linearization function of meta system for state estimation
    %
    %    varargout = linearize_joint(InferenceDS, state, V, N, U1, U2, varargin)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         N               : (c-vector) meta system observation noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %         U2              : (c-vector) meta system exogenous input 2
    %         varargin        : (strings) linearization terms wanted, e.g. 'A','B','G',....
    %    OUTPUT
    %         varargout       : (matrices) linearization terms corresponding with varargin strings

    model = InferenceDS.model;                                                     % copy existing model

    dimX  = model.statedim;
%    dimO  = model.obsdim;
%    dimN  = model.Ndim;
    dimV  = model.Vdim;

    paramIdxVec = InferenceDS.paramParamIdxVec;
    dimW  = length(paramIdxVec);
    ext_state  = state(1:dimX);
    ext_params = state(dimX+1:dimX+dimW);
    ext_V      = V(1:dimV);
%    param_V    = V(dimV+1:dimV+dimW);

    % Setup temporary model to use for linearization purposes
    model = model.setparams( model, ext_params, paramIdxVec);       % set parameters acording to state variable


    nop = length(varargin);

    for k=1:nop,

        switch varargin{k}

         case 'A'
            A = eye(InferenceDS.statedim);
            A(1:dimX,1:dimX) = model.linearize( model, ext_state, ext_V, N, U1, U2, 'A');
            A(1:dimX,dimX+1:dimX+dimW) = model.linearize( model, ext_state, ext_V, N, U1, U2, 'JFW', paramIdxVec);
            varargout{k} = A;

        case 'B'
            B = zeros(InferenceDS.statedim,InferenceDS.U1dim);
            B(1:dimX,:) = model.linearize( model, ext_state, ext_V, N, U1, U2, 'B');
            varargout{k} = B;

        case 'G'
            G = zeros(InferenceDS.statedim,InferenceDS.Vdim);
            G(1:dimX,1:dimV) = model.linearize( model, ext_state, ext_V, N, U1, U2, 'G');
            G(dimX+1:dimX+dimW,dimV+1:dimV+dimW) = eye(dimW);
            varargout{k} = G;

         case 'C'
            C = zeros(InferenceDS.obsdim,InferenceDS.statedim);
            C(:,1:dimX) =  model.linearize( model, ext_state, ext_V, N, U1, U2, 'C');
            C(:,dimX+1:dimX+dimW) = model.linearize( model, ext_state, ext_V, N, U1, U2, 'JHW', paramIdxVec);
            varargout{k} = C;

         case 'D'
            varargout{k} = model.linearize( model, ext_state, ext_V, N, U1, U2, 'D');

         case 'H'
            varargout{k} = model.linearize( model, ext_state, ext_V, N, U1, U2, 'H');

        end


    end



%-------------------------------------------------------------------------------------
function innov = innovation_generic(InferenceDS, obs, observ)

    %  INNOVATION_GENERIC  Calculates the innovation signal (difference) between the
    %   output of HFUN, i.e. OBSERV (the predicted system observation) and an actual
    %   'real world' observation OBS.
    %
    %   innov = innovation_generic(InferenceDS, obs, observ)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         obs             : (c-vector)  real-world observation vector
    %         observ          : (c-vector)  meta system generated observation vector
    %    OUTPUT
    %         inov            : (c-vector) innovation sequence

    innov = obs - observ;


%-------------------------------------------------------------------------------------
function varargout = linearize_generic(InferenceDS, state, V, N, U1, U2, varargin)

    %  LINEARIZE_GENERIC  Generic (perturbation based) linearization function of meta system
    %
    %    varargout = linearize_generic(InferenceDS, state, V, N, U1, U2, varargin)
    %
    %    INPUT
    %         InferenceDS     : (InferenceDS) Inference data structure
    %         state           : (c-vector) meta system state vector
    %         V               : (c-vector) meta system process noise vector
    %         N               : (c-vector) meta system observation noise vector
    %         U1              : (c-vector) meta system exogenous input 1
    %         U2              : (c-vector) meta system exogenous input 2
    %         varargin        : (strings) linearization terms wanted, e.g. 'A','B','G',....
    %    OUTPUT
    %         varargout       : (matrices) linearization terms corresponding with varargin strings

  nop = length(varargin);

  epsilon = 1e-8;                    % perturbation step size

  for k=1:nop,

      switch varargin{k}

      case 'A'
      %%%========================================================
      %%%             Calculate A = dffun/dstate
      %%%========================================================
      A  = zeros(InferenceDS.statedim);
      S = cvecrep(state, InferenceDS.statedim+1);
      S(:,2:end) = S(:,2:end)+diag(cvecrep(epsilon,InferenceDS.statedim));
      F = InferenceDS.ffun( InferenceDS, S, cvecrep(V,InferenceDS.statedim+1), cvecrep(U1,InferenceDS.statedim+1));
      A = (F(:,2:end)-cvecrep(F(:,1),InferenceDS.statedim))./epsilon;
      varargout{k} = A;


      case 'B'
      %%%========================================================
      %%%             Calculate B = dffun/dU1
      %%%========================================================
      %%%========================================================
      B = zeros(InferenceDS.statedim,InferenceDS.U1dim);
      U = cvecrep(U1, InferenceDS.U1dim+1);
      U(:,2:end) = U(:,2:end)+diag(cvecrep(epsilon,InferenceDS.U1dim));
      F = InferenceDS.ffun( InferenceDS, cvecrep(state,InferenceDS.U1dim+1, cvecrep(V,InferenceDS.U1dim+1), U));
      B = (F(:,2:end)-cvecrep(F(:,1),InferenceDS.U1dim))./epsilon;
      varargout{k} = B;


      case 'G'
      %%%========================================================
      %%%             Calculate G = dffun/dv
      %%%========================================================
      G = zeros(InferenceDS.statedim,InferenceDS.Vdim);
      VV = cvecrep(V, InferenceDS.Vdim+1);
      VV(:,2:end) = VV(:,2:end)+diag(cvecrep(epsilon,InferenceDS.Vdim));
      F = InferenceDS.ffun( InferenceDS, cvecrep(state,InferenceDS.Vdim+1), VV, cvecrep(U1,InferenceDS.Vdim+1));
      G = (F(:,2:end)-cvecrep(F(:,1),InferenceDS.Vdim))./epsilon;
      varargout{k} = G;


      case 'C'
      %%%========================================================
      %%%             Calculate C = dhfun/dx
      %%%========================================================
      C = zeros(InferenceDS.obsdim,InferenceDS.statedim);
      S = cvecrep(state, InferenceDS.statedim+1);
      S(:,2:end) = S(:,2:end)+diag(cvecrep(epsilon,InferenceDS.statedim));
      F = InferenceDS.hfun( InferenceDS, S, cvecrep(N,InferenceDS.statedim+1), cvecrep(U2,InferenceDS.statedim+1));
      C = (F(:,2:end)-cvecrep(F(:,1),InferenceDS.statedim))./epsilon;
      varargout{k} = C;


      case 'D'
      %%%========================================================
      %%%             Calculate D = dhfun/dU2
      %%%========================================================
      %%%========================================================
      D = zeros(InferenceDS.obsdim,InferenceDS.U2dim);
      U = cvecrep(U2, InferenceDS.U2dim+1);
      U(:,2:end) = U(:,2:end)+diag(cvecrep(epsilon,InferenceDS.U2dim));
      F = InferenceDS.hfun( InferenceDS, cvecrep(state,InferenceDS.U2dim+1, cvecrep(N,InferenceDS.U2dim+1), U));
      D = (F(:,2:end)-cvecrep(F(:,1),InferenceDS.U2dim))./epsilon;
      varargout{k} = B;

      case 'H'
      %%%========================================================
      %%%             Calculate H = dhfun/dn
      %%%========================================================
      H = zeros(InferenceDS.obsdim,InferenceDS.Ndim);
      NN = cvecrep(N, InferenceDS.Ndim+1);
      NN(:,2:end) = NN(:,2:end)+diag(cvecrep(epsilon,InferenceDS.Ndim));
      F = InferenceDS.hfun( InferenceDS, cvecrep(state,InferenceDS.Ndim+1), NN, cvecrep(U2,InferenceDS.Ndim+1));
      H = (F(:,2:end)-cvecrep(F(:,1),InferenceDS.Ndim))./epsilon;
      varargout{k} = H;

      otherwise
        error('[ linearize_generic ] Invalid linearization term requested!');

      end


  end

