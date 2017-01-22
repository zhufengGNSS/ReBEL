function [X,comp] = gmmsample(gmmDS, N)

% GMMSAMPLE  Draw N samples from the Gaussian mixture model (GMM) described by the
%            GMM data structure 'gmmDS'.
%
%   [X,comp] = gmmsample(gmmDS, N)
%
%   INPUT
%          gmmDS         Gaussian mixture model data structure with the following fields
%            .cov_type   covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'    [string]
%            .dim        data dimension  [scalar]
%            .M          number of Gaussian component densities  [scalar]
%            .weights    mixing priors (component weights) [1-by-M matrix]
%            .mu         M Gaussian component means (columns of matrix) [dim-by-M matrix]
%            .cov        covariance matrices of Gaussian components (must comply with .cov_type)
%                        [dim-by-dim-by-N matrix]
%          N             number of samples to generate [scalar]
%   OUTPUT
%          X             buffer of N samples drawn from the GMM  [dim-by-N matrix]
%          comp          component index of samples [1-by-N vector]
%
%   See also
%     GMMEVAL
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

dim   = gmmDS.dim;                         % random vector dimension
Ncomp = gmmDS.M;                           % number of component densities
w     = gmmDS.weights(:);                  % prior mixing probabilities
mu    = gmmDS.mu;                          % component means
cov   = gmmDS.cov;                         % component covariance matrices

u    = rand(1,N);

[Nc,comp] = histc(u, cumsum([0; w]));      % draw component indices according to prior
                                           % probabilities specified in gmmDS.weights
                                           % Nc = number of samples in each component
                                           % comp = index vector

X = zeros(dim, N);

% Sample each component according to the prior probabilities
switch gmmDS.cov_type

   %----------------------------------------------------------------------
   case {'full','diag'}

   for k=1:Ncomp,
      idx = find(comp==k);
      X(:,idx) = chol(cov(:,:,k))' * randn(dim,Nc(k));
   end

   %----------------------------------------------------------------------
   case {'sqrt','sqrt-diag'}

      for k=1:Ncomp,
          idx = find(comp==k);
          X(:,idx) = cov(:,:,k) * randn(dim,Nc(k));
      end

  otherwise
    error(' [ gmmsample ] Unknown covariance type!');

end

X = X + mu(:,comp);     % add in means
