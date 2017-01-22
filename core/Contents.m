% ReBEL : Recursive Bayesian Estimation Library  - Toolkit
% Version 0.2
%
% ---CORE ROUTINES---
%
% ReBEL Inference System Routines
%    consistent        -  Check ReBEL data structures for consistency.
%    convgausns        -  Convert a Gaussian noise source from one cov_type to another.
%    fixinfds          -  Make a user defined InferenceDS data structure compliant.
%    geninfds          -  Generate a InferenceDS data structure from a user specified
%                         general state space model (GSSM file).
%    gennoiseds        -  Generate a noise source data structure.
%    gensysnoiseds     -  Generate inference system noise sources (i.e. process and
%                         observation noise sources).
%    gssm              -  TEMPLATE : General state space model template. Copy and adapt
%                         for your own use.
%
% Inference Algorithms
%    cdkf              -  Central Difference Kalman Filter (SPKF family).
%    ekf               -  Extended Kalman Filter
%    gspf              -  Gaussian Sum Particle Filter
%    gmsppf            -  Gaussian Mixture Sigma-Point Particle Filter
%    kf                -  Kalman Filter (standard linear version)
%    pf                -  Generic Particle Filter (a.k.a Bootstrap or CONDENSATION)
%    sppf              -  Sigma-Point Particle Filter (Sigma-Point Filter family)
%    srcdkf            -  Square-Root Central Difference Kalman Filter (SPKF family)
%    srukf             -  Square-Root Unscented Kalman Filter (SPKF family)
%    ukf               -  Unscented Kalman Filter (SPKF family)
%
% Neural Neworks
%    mlpff             -  Feed forward a ReBEL MLP (multi-layer perceptron) neural
%                         network.
%    mlpindexgen       -  Generate 'fast unpacking' index vectors for a ReBEL MLP neural
%                         network.
%    mlpjacobian       -  Calculate neural network derivatives.
%    mlppack           -  Pack a ReBEL MLP neural network parameters (weights and biases)
%                         into a single vector.
%    mlpunpack         -  Unpack a ReBEL MLP neural network parameter vector into seperate
%                         weight and bias matrices.
%    mlpweightinit     -  Initialize the parameters of a ReBEL MLP neural network.
%
% Other Models
%    gauseval          -  Calculate the probability ( likelihood p(x|M) ) of a dataset
%                         given a multivariate Gaussian density.
%    gaussamp          -  Sample from a multivariate Gaussian density.
%    gmmfit            -  Fit/train a Gaussian mixture model (GMM) to data using EM
%    gmminitialize     -  Initiliaze a GMM (used internally by gmmfit)
%    gmmsample         -  Sample efficiently from a GMM
%    gmmprobability    -  Calculate all probabilities relating a dataset to a GMM
%                         (i.e. likelihoods, priors, evidence & posterior)
%
% Miscellaneous
%    addangle          -  Add two angles MOD 2pi radians.
%    addrelpath        -  Add and expand a relative path to the current MATLABPATH
%    checkdups         -  Check a vector for duplicate entries.
%    checkstructfields -  Check if structure has a list of specified fields.
%    cvecrep           -  Column vector replicate.
%    datamat           -  Create a datamatrix from a vector of data
%    remrelpath        -  Remove a relative path from the current MATLABPATH
%    residualresample  -  Residual resampling needed by SIR algorithms
%    rvecrep           -  Row vector replicate
%    stringmatch       -  Match one string to a cell array of others.
%    subangle          -  Subtract two angles MOD 2pi radians.
%