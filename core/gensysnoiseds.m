function [pNoise, oNoise, InferenceDS] = gensysnoiseds(InferenceDS, estimatorType, pNoiseAdaptMethod, pNoiseAdaptParams, ...
                                          oNoiseAdaptMethod, oNoiseAdaptParams)

% GENSYSNOISEDS  Generate process and observation noise data structures for a given InferenceDS data structure
%                and algorithm type. All ReBEL estimation algorithms take an inference data structure (InferenceDS),
%                as well as two system noise data structures (process noise and observation noise) as arguments.
%
%   [pNoise, oNoise] = gensysnoiseds(InferenceDS, estimatorType, pNoiseAdaptMethod, pNoiseAdaptParams, oNoiseAdaptMethod, oNoiseAdaptParams))
%
%   INPUT
%          InferenceDS         (InferenceDS) Inference data structure generated from a GSSM file by 'geninfds'
%          estimatorType       (string) type of estimator to be used (i.e. 'kf', 'ukf', 'ekf', 'pf', etc.)
%          pNoiseAdaptMethod  <<optional>> (string) Process noise covariance adaptation method :
%                                      'anneal'        : annealing
%                                      'lambda-decay'  : RLS like lambda decay
%                                      'robbins-monro' : Robbins-Monro stochastic approximation
%                               If this field is set, then pNoiseAdaptParams must also be set.
%          pNoiseAdaptParams  <<optional>> (vector) noise adaptation parameters. Depend on pNoiseAdaptMethod
%                                 if 'anneal'        : [annealing_factor minimum_allowed_variance]
%                                 if 'lambda-decay'  : [lambda_factor minimum_allowed_variance]
%                                 if 'robbins-monro' : [1/nu_initial 1/nu_final]
%          oNoiseAdaptMethod  <<optional>> Observation noise covariance adaptation method : same as above
%                                          except the only allowed method is 'robbins-monro'
%          oNoiseAdaptParams  <<optional>> Same as above for process noise
%
%   OUTPUT
%          pNoise              (NoiseDS) process noise data structure
%          oNoise              (NoiseDS) observation noise data structure
%          InferenceDS         (InferenceDS) updated inference data structure
%
%     See also
%     GENINFDS, GENNOISEDS
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

%=== ERROR CHECKING ==========================================================================

if ((nargin < 2) | rem(nargin,2))
    error(' [ gensysnoiseds ] Not enough input parameters.');
end

if (nargout ~= 3)
    error(' [ gensysnoiseds ] Not enough output arguments.');
end


error(consistent(InferenceDS,'InferenceDS'));         %-- check for consistency of InferenceDS data structure

InferenceDS.esttype = estimatorType;                  % store estimator type


%=== INFERENCE TYPE SPECIFIC STRUCTURE ==================================================================

switch (InferenceDS.inftype)

%----------------------------------------- STATE ESTIMATION ---------------------------------------------
case 'state'

    %--- Generate/convert or copy noise sources from GSSM data structure

    pNoise = InferenceDS.model.pNoise;                            % process noise data structure
    oNoise = InferenceDS.model.oNoise;                            % observation noise data structure


    %--- KALMAN FILTER FAMILY : Checks and conversion

    if stringmatch(estimatorType, {'kf','ekf','ukf','cdkf','srukf','srcdkf'})

        % If default noise source is not Guassian, define a Gaussian noise source with the same dimension, mean and covariance
        % if available
        if ~stringmatch(pNoise.ns_type, {'gaussian','combo-gaussian'})
            Arg.type = 'gaussian';         % standard Gaussian noise source
            Arg.cov_type = 'full';          % with full covariance matrix
            Arg.dim = pNoise.dim;          % process noise dimension
            if isfield(pNoise,'mu')
              Arg.mu = pNoise.mu;
            else
              %warning(' [ gensysnoiseds ] Process noise data structure does not have a defined mean vector. Default assigned.');
              Arg.mu = zeros(Arg.dim,1);     % default : zero mean
            end
            if isfield(pNoise,'cov')
              Arg.cov = pNoise.cov;
            else
              %warning(' [ gensysnoiseds ] Process noise data structure does not have a defined covariance matrix. Default assigned.');
              Arg.cov  = eye(Arg.dim);         % default : covariance
            end
            pNoise = gennoiseds(Arg);      % generate process noise data structure
        end
        if ~stringmatch(oNoise.ns_type, {'gaussian','combo-gaussian'})
            Arg.type = 'gaussian';         % standard Gaussian noise source
            Arg.cov_type = 'full';          % with full covariance matrix
            Arg.dim = oNoise.dim;          % process noise dimension
            if isfield(oNoise,'mu')
              Arg.mu = oNoise.mu;
            else
              %warning(' [ gensysnoiseds ] Observation noise data structure does not have a defined mean vector. Default assigned.');
              Arg.mu = zeros(Arg.dim,1);     % default : zero mean
            end
            if isfield(oNoise,'cov')
              Arg.cov = oNoise.cov;
            else
              %warning(' [ gensysnoiseds ] Observation noise data structure does not have a defined covariance matrix. Default assigned.');
              Arg.cov  = eye(Arg.dim);         % default : covariance
            end
            oNoise = gennoiseds(Arg);      % generate observation noise data structure
        end

        if stringmatch(estimatorType, {'srukf','srcdkf'})   % Check for square root Kalman algorithms
            %-- process noise
            switch (pNoise.cov_type)                        % Determine cov_type of Gaussian noise source
            case 'diag'
                pNoise = convgausns(pNoise,'sqrt-diag');
                %warning(' [ gensysnoiseds ] Converting process noise source covariance type to ''sqrt-diag''.');
            case 'full'
                pNoise = convgausns(pNoise,'sqrt');
                %warning(' [ gensysnoiseds ] Converting process noise source covariance type to ''sqrt''.');
            end
            %-- observation noise
            switch (oNoise.cov_type)                        % Determine cov_type of Gaussian noise source
            case 'diag'
                oNoise = convgausns(oNoise,'sqrt-diag');
                %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt-diag''.');
            case 'full'
                oNoise = convgausns(oNoise,'sqrt');
                %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt''.');
            end
        end

        InferenceDS.InovUpdateMaskIdxVec = []; % Innovation update mask index vector ... Indicates which components
                                               % of the innovation vector should be ignored when calculating a Kalman
                                               % state update

    end


    %--------------------------------------------------------------------------------------------
    %--- PARTICLE FILTER FAMILY : Checks and conversion


    %----------------------------------------------------------------------------------
    if stringmatch(estimatorType, 'gspf')   % 'Gaussian Sum Particle Filter'

        % If process noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
        % if available
        if ~stringmatch(pNoise.ns_type, 'gmm')

            Arg.type = 'gmm';              % GMM noise source
            Arg.cov_type = 'sqrt';         % GSPF use square-root covariance matrices
            Arg.dim = pNoise.dim;          % process noise dimension
            Arg.M = 1;                     % single component
            Arg.weights = [1];             % component weight
            if isfield(pNoise,'mu')
              Arg.mu = pNoise.mu;
            else
              %warning(' [ gensysnoiseds ] Process noise data structure does not have a defined mean vector. Default assigned');
              Arg.mu = zeros(Arg.dim,1);     % default : zero mean
            end
            Arg.cov = repmat(zeros(Arg.dim),[1 1 1]); % default covariance buffer

            if isfield(pNoise,'cov')
              if isfield(pNoise,'cov_type') cov_type = pNoise.cov_type; else cov_type='full'; end
              switch cov_type
              case {'sqrt','sqrt-diag'}
                  Arg.cov(:,:,1) = pNoise.cov;
              case {'full','diag'}
                  Arg.cov(:,:,1) = chol(pNoise.cov)';
              otherwise
                  error('[ gensysnoiseds::gspf ] Unknown process noise covariance type.');
              end
            else
              warning(' [ gensysnoiseds::gspf ] Process noise data structure does not have a defined covariance matrix. Default assigned.');
              Arg.cov(:,:,1) = repmat(eye(Arg.dim),[1 1 1]);  % default : covariance  (Cholesky factor)
            end
            pNoise = gennoiseds(Arg);      % generate process noise data structure

        else

        % Make sure the GMM component densities is of cov_type 'sqrt'

            %-- process noise
            switch (pNoise.cov_type)                        % Determine cov_type of Gaussian noise source
            case 'diag'
                pNoise = convgausns(pNoise,'sqrt-diag');
                %warning(' [ gensysnoiseds ] Converting process noise source covariance type to ''sqrt-diag''.');
            case 'full'
                pNoise = convgausns(pNoise,'sqrt');
                %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt''.');
            end

        end

    end


    %----------------------------------------------------------------------------------
    if stringmatch(estimatorType, 'gmsppf')   % 'Gaussian Mixture Sigma-Point Particle Filter'

        % If process noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
        % if available
        if ~stringmatch(pNoise.ns_type, 'gmm')
            Arg.type = 'gmm';              % GMM noise source
            Arg.cov_type = 'sqrt';         % GSPF use square-root covariance matrices
            Arg.dim = pNoise.dim;          % process noise dimension
            Arg.M = 1;                     % single component
            Arg.weights = [1];             % component weight
            if isfield(pNoise,'mu')
              Arg.mu = pNoise.mu;
            else
              warning(' [ gensysnoiseds ] Process noise data structure does not have a defined mean vector. Default assigned');
              Arg.mu = zeros(Arg.dim,1);     % default : zero mean
            end
            Arg.cov = repmat(zeros(Arg.dim),[1 1 1]); % default covariance buffer
            if isfield(pNoise,'cov')
              if isfield(pNoise,'cov_type') cov_type = pNoise.cov_type; else cov_type='full'; end
              switch cov_type
                case {'sqrt','sqrt-diag'}
                    Arg.cov(:,:,1) = pNoise.cov;
                case {'full','diag'}
                    Arg.cov(:,:,1) = chol(pNoise.cov)';
                otherwise
                    error(' [ gensysnoiseds::gmsppf ] Unknown process noise covariance type.');
              end

            else
              warning(' [ gensysnoiseds::gmsppf ] Process noise data structure does not have a defined covariance matrix. Default assigned.');
              Arg.cov(:,:,1) = repmat(eye(Arg.dim),[1 1 1]);  % default : covariance  (Cholesky factor)
            end
            pNoise = gennoiseds(Arg);      % generate process noise data structure
        else
            % Make sure the GMM component densities is of cov_type 'sqrt'
            %-- process noise
            switch (pNoise.cov_type)                        % Determine cov_type of Gaussian noise source
            case 'diag'
                pNoise = convgausns(pNoise,'sqrt-diag');
                %warning(' [ gensysnoiseds ] Converting process noise source covariance type to ''sqrt-diag''.');
            case 'full'
                pNoise = convgausns(pNoise,'sqrt');
                %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt''.');
            end
        end

        % If observation noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
        % if available
        if ~stringmatch(oNoise.ns_type, 'gmm')
            Arg.type = 'gmm';              % GMM noise source
            Arg.cov_type = 'sqrt';         % GSPF use square-root covariance matrices
            Arg.dim = oNoise.dim;          % observation noise dimension
            Arg.M = 1;                     % single component
            Arg.weights = [1];             % component weight
            if isfield(oNoise,'mu')
              Arg.mu = oNoise.mu;
            else
              warning(' [ gensysnoiseds::gmsppf ] Observation noise data structure does not have a defined mean vector. Default assigned');
              Arg.mu = zeros(Arg.dim,1);     % default : zero mean
            end
            Arg.cov = repmat(zeros(Arg.dim),[1 1 1]); % default covariance buffer
            if isfield(oNoise,'cov')
              if isfield(oNoise,'cov_type') cov_type = oNoise.cov_type; else cov_type='full'; end
              switch cov_type
              case {'sqrt','sqrt-diag'}
                  Arg.cov(:,:,1) = oNoise.cov;
              case {'full','diag'}
                  Arg.cov(:,:,1) = chol(oNoise.cov)';
              otherwise
                  error(' [ gensysnoiseds::gmsppf ] Unknown observation noise covariance type.');
              end
            else
              warning(' [ gensysnoiseds::gmsppf ] Observation noise data structure does not have a defined covariance matrix. Default assigned.');
              Arg.cov(:,:,1) = repmat(eye(Arg.dim),[1 1 1]);  % default : covariance  (Cholesky factor)
            end
            oNoise = gennoiseds(Arg);      % generate process noise data structure
        else
            % Make sure the GMM component densities is of cov_type 'sqrt'
            %-- process noise
            switch (oNoise.cov_type)                        % Determine cov_type of Gaussian noise source
            case 'diag'
                oNoise = convgausns(oNoise,'sqrt-diag');
                %warning(' [ gensysnoiseds ] Converting process noise source covariance type to ''sqrt-diag''.');
            case 'full'
                oNoise = convgausns(oNoise,'sqrt');
                %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt''.');
            end
        end

    end


    %--- Setup noise source tags

    pNoise.tag = 'state';                   % tag this as a state variable noise sources
    oNoise.tag = 'obs';                     % tag this as a observation variable noise source



%----------------------------------------- PARAMETER ESTIMATION ---------------------------------------------
case 'parameter'

    %--- Generate default process noise source

    switch estimatorType
    case {'gspf','gmsppf'}
        Arg.type = 'gmm';                                    % standard Gaussian noise source
        Arg.cov_type = 'sqrt';                               % with full covariance matrix
        Arg.tag = 'param';                                   % this noise source operates on parameters
        Arg.dim = InferenceDS.Vdim;                          % process noise dimension
        Arg.mu = zeros(Arg.dim,1);                           % default : zero mean
        Arg.M  = 1;                                          % single component GMM
        Arg.weights = [1];                                   % component weight
        Arg.cov  = repmat(eye(Arg.dim),[1 1 1]);             % default : unity covariance
    otherwise
        Arg.type = 'gaussian';                               % standard Gaussian noise source
        Arg.cov_type = 'full';                               % with full covariance matrix
        Arg.tag = 'param';                                   % this noise source operates on parameters
        Arg.dim = InferenceDS.Vdim;                          % process noise dimension
        Arg.mu = zeros(Arg.dim,1);                           % default : zero mean
        Arg.cov  = eye(Arg.dim);                             % default : unity covariance
    end

    pNoise = gennoiseds(Arg);                                 % generate default process noise source


    %--- Generate default observation noise source


    %--- KALMAN FILTER FAMILY : Checks and conversion

    if stringmatch(estimatorType, {'kf','ekf','ukf','cdkf','srukf','srcdkf'})

        if stringmatch(InferenceDS.paramFunSelect, {'both','both-p','ffun'})
            % If default noise source is not Guassian, define a Gaussian noise source with the same dimension
            if ~stringmatch(InferenceDS.model.pNoise.ns_type, {'gaussian','combo-gaussian'})
                Arg.type = 'gaussian';         % standard Gaussian noise source
                Arg.cov_type = 'full';          % with full covariance matrix
                Arg.dim = InferenceDS.model.pNoise.dim;          % process noise dimension
                if isfield(InferenceDS.model.pNoise,'mu')
                   Arg.mu = InferenceDS.model.pNoise.mu;
                else
                   Arg.mu = zeros(Arg.dim,1);     % default : zero mean
                   %warning(' [ gensysnoiseds ] Process noise data structure does not have a defined mean vector. Default assigned.');
                end
                if isfield(InferenceDS.model.pNoise,'cov')
                   Arg.cov = InferenceDS.model.pNoise.cov;
                else
                   Arg.cov = eye(Arg.dim);         % default : covariance
                   %warning(' [ gensysnoiseds ] Process noise data structure does not have a defined covariance matrix. Default assigned.');
                end
                obs_pNoise = gennoiseds(Arg);      % generate process noise data structure
            else
                obs_pNoise = InferenceDS.model.pNoise; % or copy the original pnoise if its already gaussian
            end
        end
        if stringmatch(InferenceDS.paramFunSelect, {'both','both-p','hfun'})
            % If default noise source is not Guassian, define a Gaussian noise source with the same dimension
            if ~stringmatch(InferenceDS.model.oNoise.ns_type, {'gaussian','combo-gaussian'})
                Arg.type = 'gaussian';         % standard Gaussian noise source
                Arg.cov_type = 'full';          % with full covariance matrix
                Arg.dim = InferenceDS.model.oNoise.dim;          % process noise dimension
                if isfield(InferenceDS.model.oNoise,'mu')
                    Arg.mu = InferenceDS.model.oNoise.mu;
                else
                    Arg.mu = zeros(Arg.dim,1);     % default : zero mean
                    %warning(' [ gensysnoiseds ] Observation noise data structure does not have a defined mean vector. Default assigned.');
                end
                if isfield(InferenceDS.model.oNoise,'cov')
                    Arg.cov = InferenceDS.model.oNoise.cov;
                else
                    Arg.cov = eye(Arg.dim);         % default : covariance
                    %warning(' [ gensysnoiseds ] Observation noise data structure does not have a defined covariance matrix. Default assigned.');
                end
                obs_oNoise = gennoiseds(Arg);      % generate observation noise data structure
            else
                obs_oNoise = InferenceDS.model.oNoise; % or copy the original onoise if its already gaussian
            end
        end

    else

        if stringmatch(InferenceDS.paramFunSelect, {'both','both-p','ffun'})
            obs_pNoise = InferenceDS.model.pNoise; % or copy the original pnoise if its already gaussian
        end
        if stringmatch(InferenceDS.paramFunSelect, {'both','both-p','hfun'})
            obs_oNoise = InferenceDS.model.oNoise; % or copy the original onoise if its already gaussian
        end

    end


    %--- NOW BUILD OBSERVATION NOISE SOURCE

    switch InferenceDS.paramFunSelect
    case {'both','both-p'}
        clear Arg;
        Arg.tag = 'obs';                        % ID tag : this noise source operates on observations
        Arg.dim = InferenceDS.Ndim;             % set noise source dimension
        Arg.type = 'combo';                     % Combination noise source
        Arg.noiseSources = {obs_pNoise , obs_oNoise};   % construct noise source cell array
        oNoise = gennoiseds(Arg);               % generate observation noise source
     case 'ffun'
        oNoise = obs_pNoise;                    % copy process noise source from GSSM
        oNoise.tag = 'obs';                     %  ID tag : this noise source operates on observations
     case 'hfun'
        oNoise = obs_oNoise;                    % copy observation noise source from GSSM
        oNoise.tag = 'obs';                     %  ID tag : this noise source operates on observations
     otherwise
        error([' [ gensysnoiseds::parameter ] Unknown paramFunSelect value ''' InferenceDS.paramFunSelect '''']);
    end


    %--- KALMAN FILTER FAMILY : Checks and conversion (again!)

    if stringmatch(estimatorType, {'srukf','srcdkf'})           % Check for square root Kalman algorithms
        %-- process noise
        pNoise = convgausns(pNoise,'sqrt');
        %-- observation noise
        switch (oNoise.cov_type)                                % Determine cov_type of Gaussian noise source
        case 'diag'
            oNoise = convgausns(oNoise,'sqrt-diag');
            %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt-diag''.');
        case 'full'
            oNoise = convgausns(oNoise,'sqrt');
            %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt''.');
        end
    end




%----------------------------------------- JOINT ESTIMATION ---------------------------------------------
case 'joint'

    %--- Generate/convert or copy noise sources from GSSM data structure

    param_pNoise_Arg.type = 'gaussian';              % standard Gaussian noise source
    param_pNoise_Arg.cov_type = 'full';              % with full covariance matrix
    param_pNoise_Arg.tag = 'param';                  % this noise source operates on parameters
    param_pNoise_Arg.dim = length(InferenceDS.paramParamIdxVec); % noise dimension (length of parameter vector)
    param_pNoise_Arg.mu  = zeros(param_pNoise_Arg.dim,1);         % default : zero mean
    param_pNoise_Arg.cov = eye(param_pNoise_Arg.dim);             % default : unity covariance

    state_pNoise = InferenceDS.model.pNoise;        % Copy GSSM process noise source for state part of state vector
    state_pNoise.tag = 'state';

    oNoise = InferenceDS.model.oNoise;              % The observation noise source is the same as that of the
                                                    % underlying model

    %--- KALMAN FILTER FAMILY : Checks and conversion

    if stringmatch(estimatorType, {'kf','ekf','ukf','cdkf','srukf','srcdkf'})

        % If default noise source is not Guassian, define a Gaussian noise source with the same dimension
        if ~stringmatch(state_pNoise.ns_type, {'gaussian','combo-gaussian'})
            clear Arg;
            Arg.type = 'gaussian';          % standard Gaussian noise source
            Arg.cov_type = 'full';           % with full covariance matrix
            Arg.dim = state_pNoise.dim;     % process noise dimension
            if isfield(state_pNoise,'mu')
              Arg.mu = state_pNoise.mu;
            else
              Arg.mu = zeros(Arg.dim,1);     % default : zero mean
              %warning(' [ gensysnoiseds ] Process noise data structure does not have a defined mean vector. Default assigned.');
            end
            if isfield(state_pNoise,'cov')
              Arg.cov = state_pNoise.cov;
            else
              Arg.cov  = eye(Arg.dim);         % default : covariance
              %warning(' [ gensysnoiseds ] Process noise data structure does not have a defined covariance matrix. Default assigned.');
            end
            state_pNoise = gennoiseds(Arg); % generate process noise data structure
        else
            param_pNoise_Arg.cov_type = state_pNoise.cov_type;
        end

        if ~stringmatch(oNoise.ns_type, {'gaussian','combo-gaussian'})
            Arg.type = 'gaussian';         % standard Gaussian noise source
            Arg.cov_type = 'full';          % with full covariance matrix
            Arg.dim = oNoise.dim;          % process noise dimension
            if isfield(oNoise,'mu')
              Arg.mu = oNoise.mu;
            else
              Arg.mu = zeros(Arg.dim,1);     % default : zero mean
              %warning(' [ gensysnoiseds ] Observation noise data structure does not have a defined mean vector. Default assigned.');
            end
            if isfield(oNoise,'cov')
              Arg.cov = oNoise.cov;
            else
              Arg.cov  = eye(Arg.dim);         % default : covariance
              %warning(' [ gensysnoiseds ] Observation noise data structure does not have a defined covariance matrix. Default assigned.');
            end
            oNoise = gennoiseds(Arg);      % generate observation noise data structure
        end

    end


    %--- PARTICLE FILTER FAMILY : Checks and conversions

    if stringmatch(estimatorType, {'gspf','gmsppf'})

        error(' [ gensysnoiseds ] Joint estimation is not yet supported in the meta/abstract level for algorithms : GSPF');

        %param_pNoise_Arg.type = 'gmm';              %
        %param_pNoise_Arg.cov_type = 'sqrt';         % with full covariance matrix
        %param_pNoise_Arg.tag = 'param';                  % this noise source operates on parameters
        %param_pNoise_Arg.dim = length(InferenceDS.paramParamIdxVec); % noise dimension (length of parameter vector)
        %param_pNoise_Arg.mu  = zeros(param_pNoise_Arg.dim,1);         % default : zero mean
        %param_pNoise_Arg.cov = eye(param_pNoise_Arg.dim);             % default : unity covariance

        %state_pNoise = InferenceDS.model.pNoise;        % Copy GSSM process noise source for state part of state vector
        %state_pNoise.tag = 'state';

        %oNoise = InferenceDS.model.oNoise;              % The observation noise source is the same as that of the
                                                    % underlying model


    end


    param_pNoise = gennoiseds(param_pNoise_Arg);   % Generate default process noise source for parameter part of state vector


    clear Arg;
    Arg.tag = 'state/param';                    % set descriptive tag
    Arg.dim = InferenceDS.Vdim;                 % set noise source dimension
    Arg.type = 'combo';                         % Combination noise source
    Arg.noiseSources = {state_pNoise , param_pNoise};   % construct noise source cell array

    pNoise = gennoiseds(Arg);                   % Generate process noise source as a combination of the original model
                                                % process noise source and the above generated artifical process noise
                                                % source operating on the parameter values in the state vector


    %--- KALMAN FILTER FAMILY : Checks and conversion   (again!)

    if stringmatch(estimatorType, {'srukf','srcdkf'})           % Check for square root algorithms
        %-- process noise
        switch (pNoise.cov_type)                                % Determine cov_type of Gaussian noise source
        case 'diag'
            pNoise = convgausns(pNoise,'sqrt-diag');
            %warning(' [ gensysnoiseds ] Converting process noise source covariance type to ''sqrt-diag''.');
        case 'full'
            pNoise = convgausns(pNoise,'sqrt');
            %warning(' [ gensysnoiseds ] Converting process noise source covariance type to ''sqrt''.');
        end
        %-- observation noise
        switch (oNoise.cov_type)                                % Determine cov_type of Gaussian noise source
        case 'diag'
            oNoise = convgausns(oNoise,'sqrt-diag');
            %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt-diag''.');
        case 'full'
            oNoise = convgausns(oNoise,'sqrt');
            %warning(' [ gensysnoiseds ] Converting observation noise source covariance type to ''sqrt''.');
        end
    end


    %--- Setup noise source tags

    pNoise.tag = 'state/param';             % tag this as a state variable noise sources
    oNoise.tag = 'obs';                     % tag this as a observation variable noise source



%--------------------------------------------------------------------------------------------------------
otherwise

  error([' [ gensysnoiseds ] Unknown inference type ''' InferenceDS.inftype ''' in InferenceDS.type ']);

end

InferenceDS.pNoiseAdaptMethod = [];
InferenceDS.oNoiseAdaptMethod = [];
pNoise.adaptMethod = [];
oNoise.adaptMethod = [];

if (nargin >= 4)
  pNoise.adaptMethod = pNoiseAdaptMethod;
  pNoise.adaptParams = pNoiseAdaptParams;
  InferenceDS.pNoiseAdaptMethod = pNoiseAdaptMethod;
end
if (nargin == 6)
  oNoise.adaptMethod = oNoiseAdaptMethod;
  oNoise.adaptParams = oNoiseAdaptParams;
  InferenceDS.oNoiseAdaptMethod = oNoiseAdaptMethod;
end


%=== Other default parmeters
if stringmatch(estimatorType, {'ekf'})
  if ~isfield(InferenceDS,'ekfParams')
    InferenceDS.ekfParams = 1;
  end
end
