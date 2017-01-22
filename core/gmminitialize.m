function gmmDS = gmminitialize(gmmDS, X, maxI)

% GMMINITIALIZE  Initialises Gaussian mixture model (GMM) from data
%
%   INPUT
%          gmmDS         Gaussian mixture model data structure with the following fields
%            .cov_type   covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'    [string]
%            .dim        data dimension  [scalar]
%            .M          number of Gaussian component densities  [scalar]
%            .weights    mixing priors (component weights) [1-by-M matrix]
%            .mu         N Gaussian component means (columns of matrix) [dim-by-N matrix]
%            .cov        covariance matrices of Gaussian components (must comply with .cov_type)
%                        [dim-by-dim-by-N matrix]
%          X             dataset of M samples (column vectors) [dim-by-M matrix]
%          maxI          (optional) maximum number of iterations (default = 100)
%
%   OUTPUT
%          gmmDS         data initilized (updated) GMM data structure
%
%   See also
%   GMMEVAL, GMMSAMPLE
%

%   This function has been derived and modified from the 'gmminit' function in
%   the NETLAB toolkit (by Ian T Nabney and Chris Bishop). See LICENSE file
%   in the NETLAB subdirectory for the Netlab license
%   Copyright (c) Ian T Nabney (1996-2001)
%
%   The license for the derived file (this function) follows:
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

if nargin < 3
  maxI = 100;
end

[xdim, nov] = size(X);

dim = xdim;         % state dimension
M   = gmmDS.M;      % number of components

X=X';

% Arbitrary width used if variance collapses to zero: make it 'large' so
% that centre is responsible for a reasonable number of points.
GMM_WIDTH = 1.0;

  % Use kmeans algorithm to initialise the centroids from the data
  options = foptions;
  options(1) = -1;       % don't display warnings
  options(14) = maxI;    % Just use 5 iterations of k-means in initialisation
  options(5)  = 1;       % initilize centroids and their covariances from data
  [mu, options, post] = kmeans(gmmDS.mu', X, options);       % call Netlab k-means algorithm

  gmmDS.mu = mu';   % convert from Netlab format to ReBEL format

  % Set priors depending on number of points in each cluster
  cluster_sizes = max(sum(post, 1), 1);             % Make sure that no prior is zero
  gmmDS.weights = cluster_sizes/sum(cluster_sizes); % Normalise priors


fixCov = GMM_WIDTH*eye(dim);

switch gmmDS.cov_type

  case 'full'
    for j = 1:M
      % Pick out data points belonging to this centre
      c = X(find(post(:, j)),:);
      sizec = size(c,1);
      tmu = mu(j,:);
      diffs = c - tmu(ones(1,sizec),:);
      gmmDS.cov(:,:,j) = (diffs'*diffs)/sizec;
      % Add GMM_WIDTH*Identity to rank-deficient covariance matrices
      if rank(gmmDS.cov(:,:,j)) < dim
          gmmDS.cov(:,:,j) = gmmDS.cov(:,:,j) + fixCov;
      end
    end

  case 'diag'
    for j = 1:M
      % Pick out data points belonging to this centre
      c = X(find(post(:, j)),:);
      sizec = size(c,1);
      tmu = mu(j,:);
      diffs = c - tmu(ones(1,sizec),:);
      d = sum((diffs.*diffs), 1)/sizec;
      % Replace small entries by GMM_WIDTH value
      d = d + GMM_WIDTH*(d<eps);
      gmmDS.cov(:,:,j) = diag(d);
    end

  case 'sqrt'
    for j = 1:M
      % Pick out data points belonging to this centre
      c = X(find(post(:, j)),:);
      sizec = size(c,1);
      tmu = mu(j,:);
      diffs = c - tmu(ones(1,sizec),:);
      cov = (diffs'*diffs)/sizec;
      % Add GMM_WIDTH*Identity to rank-deficient covariance matrices
      if rank(cov) < gmmDS.dim
          cov = cov + fixCov;
      end
      gmmDS.cov(:,:,j) = chol(cov)';
    end

  case 'sqrt-diag'
    for j = 1:M
      % Pick out data points belonging to this centre
      c = x(find(post(:, j)),:);
      sizec = size(c,1);
      tmu = mu(j,:);
      diffs = c - tmu(ones(1,sizec),:);
      d = sum((diffs.*diffs), 1)/sizec;
      % Replace small entries by GMM_WIDTH value
      d = d + GMM_WIDTH*(d<eps);
      gmmDS.cov(:,:,j) = diag(sqrt(d));
    end


  otherwise
    error([' [ gmminitialize ] Unknown covariance type ', gmmDS.cov_type]);

end

