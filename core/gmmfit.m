function [gmmDS, leb] = gmmfit(X, M, tt, cov_type, check_cov, display, W)

% GMMFIT   Fit a Gaussian mixture model (GMM) with M components to dataset X
%          using an expectation maximization algorithm (EM).
%
%
%   [gmmDS, leb] = gmmfit(X, M, tt, cov_type,  check_cov, display, W)
%
%   INPUT
%          X            Dataset of N samples (column vectors) [dim-by-N matrix]
%          M            Number of Gaussian mixture component densities [scalar]
%                       OR a pre-initialized GMM data structure (gmmDS)
%          tt           Termination threshold 0 < tt < 1 (if % change in log likelihood
%                       falls below this value, the EM algorithm terminates.) [scalar]
%                       OR if tt is a [1-by-2 vector] the first component is the termination
%                       threshold and the second component is the maximum number of iterations
%                       allowed for the EM, i.e. tt = [tt max_iterations]
%          cov_type     Covariance type 'full','diag','sqrt','sqrt-diag' [string]
%          check_cov    (optional) Covariance check flag : If this flag is set, a covariance
%                       matrix is reset to its original value when any of its singular values
%                       are too small (less than MIN_COVAR which has the value eps). With the
%                       default value of 0 no action is taken.
%          display      (optional) Display results of training if this is set to 1 (default=0)
%          W            (optional) 1xN vector of sample weights used to do a weighted EM fit to
%                                  the data.
%
%   OUTPUT
%          gmmDS         Gaussian mixture model data structure with the following fields
%            .cov_type   covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'    [string]
%            .dim        data dimension  [scalar]
%            .M          number of Gaussian component densities  [scalar]
%            .weights    mixing priors (component weights) [1-by-M matrix]
%            .mu         N Gaussian component means (columns of matrix) [dim-by-N matrix]
%            .cov        covariance matrices of Gaussian components (must comply with .cov_type)
%                        [dim-by-dim-by-N matrix]
%          leb           log evidence buffer (sum of log evidence of data as a function of EM iteration)
%
%
%   See also
%     GMMEVAL, GMMSAMPLE, GMMINITIALIZE

%   This function has been derived and modified from the 'gmmem' function in
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
Nin = nargin;

switch Nin
case 2
  tt = [0.1 100];
  cov_type = 'full';
  check_cov = 0;
  display = 0;
case 3
  cov_type = 'full';
  check_cov = 0;
  display = 0;
case 4
  check_cov = 0;
  display = 0;
case 5
  display = 0;
case {6,7}
otherwise
  error(' [ gmmfit ] Incorrect number of input arguments.');
end

[dim, Ndata] = size(X);             % get dimensions of dataset

% Sort out the options
if (length(tt)==2)                  % number of EM iterations
  niters = tt(2);
else
  niters = 100;
end


store = 0;
if (nargout > 1)
  store = 1;    % Store the evidence values to return them
  leb = zeros(1, niters);
end

test = 0;
if tt(1) > 0.0
  test = 1;                         % Test log likelihood for termination
  leb = zeros(1, niters);
end


if isnumeric(M)                    % Initialize GMM from data if needed
  gmmDS.cov_type = cov_type;
  gmmDS.dim = dim;
  gmmDS.M = M;
  gmmDS.weights = ones(dim,M);
  gmmDS.mu = zeros(dim,M);
  gmmDS.cov = repmat(eye(dim),[1 1 M]);   % initialize GMM centroids and their covariances
  gmmDS = gmminitialize(gmmDS, X, 5);     % using Netlab's KMEANS algorithm.
else
  gmmDS = M;
  M = gmmDS.M;
end


if check_cov                        % Ensure that covariances don't collapse
  MIN_COVAR = eps;                  % Minimum singular value of covariance matrix
  MIN_COVAR_SQRT = sqrt(eps);
  init_covars = gmmDS.cov;
end


eold = -Inf;

ones_dim = ones(1,dim);
ones_Ndata = ones(Ndata,1);

% Main loop of algorithm
for n = 1:niters

  % Calculate posteriors based on old parameters
  if (nargin == 7)
      [prior, likelihood, evidence, posterior] = gmmprobability(gmmDS, X, W);
  else
      [prior, likelihood, evidence, posterior] = gmmprobability(gmmDS, X);
  end

  % Calculate error value if needed
  if (display | store | test)
    % Error value is negative log likelihood of data evidence
    e = sum(log(evidence));
    if store
      leb(n) = e;
    end
    if display
      fprintf(1, 'Cycle %4d  Evidence %11.6f\n', n, e);
    end
    if test
      if (n > 1 & abs((e - eold)/eold) < tt(1))
        leb=leb(1:n);
        return;
      else
        eold = e;
      end
    end
  end


  % Adjust the new estimates for the parameters
  new_pr = (sum(posterior, 2))';
  new_c =  (posterior * X')';

  % Now move new estimates to old parameter vectors
  gmmDS.weights = new_pr / Ndata;
  new_mu = new_c ./ new_pr(ones_dim,:);
  gmmDS.mu = new_mu;

  switch cov_type

  case 'full'
    for j = 1:M
      tmu = new_mu(:,j);
      diffs = X - tmu(:,ones_Ndata);
      tpost = sqrt(posterior(j,:));
      diffs = diffs .* tpost(ones_dim,:);
      gmmDS.cov(:,:,j) = (diffs*diffs')/new_pr(j);
    end
    if check_cov
      % Ensure that no covariance is too small
      for j = 1:M
        if min(svd(gmmDS.cov(:,:,j))) < MIN_COVAR
          gmmDS.cov(:,:,j) = init_covars(:,:,j);
        end
      end
    end

  case {'sqrt','sqrt-diag'}
    for j = 1:M
      tmu = new_mu(:,j);
      diffs = X - tmu(:,ones_Ndata);
      tpost = (1/sqrt(new_pr(j))) * sqrt(posterior(j,:));
      diffs = diffs .* tpost(ones_dim,:);
      [foo,tcov] = qr(diffs',0);
      gmmDS.cov(:,:,j) = tcov';
    end
    if check_cov
      % Ensure that no covariance is too small
      for j = 1:M
        if min(abs(diag(gmmDS.cov(:,:,j)))) < MIN_COVAR_SQRT
          gmmDS.cov(:,:,j) = init_covars(:,:,j);
        end
      end
    end


  case 'diag'
    for j = 1:M
      tmu = new_mu(:,j);
      diffs = X - tmu(:,ones_Ndata);
      tpost = posterior(j,:);
      dd = sum((diffs.*diffs).*tpost(ones_dim,:), 2)/new_pr(j);
      if check_cov
        if min(dd) < MIN_COVAR
          gmmDS.cov(:,:,j) = init_covars(:,:,j);
        else
          gmmDS.cov(:,:,j) = diag(dd);
        end
      else
        gmmDS.cov(:,:,j) = diag(dd);
      end
    end


   otherwise
      error(['Unknown covariance type ', cov_type]);
  end

end

