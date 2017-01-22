function [prior, likelihood, evidence, posterior] = gmmprobability(gmmDS, X, W)

% GMMPROBABILITY  Calculates any of the related (through Bayes rule) probabilities
%                 of a Gaussian Mixture Model (gmmDS) and a given dataset X.
%                 'prob_type' is a string indicating which of the four probability
%                 values are needed. These probabilities are:
%
%                      P(X|C) . P(C)                       likelihood . prior
%           P(C|X) = -----------------       posterior =  --------------------
%                          P(X)                                evidence
%
%            where C is the component classes (Gaussians) of the GMM and X is the data.
%
%   probability = gmmprobability(gmmDS, X, W)
%
%   INPUT
%          gmmDS         Gaussian mixture model data structure with the following fields
%            .cov_type   covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'    [string]
%            .dim        data dimension  [scalar]
%            .M          number of Gaussian component densities  [scalar]
%            .weights    mixing priors (component weights) [1-by-M vector]
%            .mu         M Gaussian component means (columns of matrix) [dim-by-M matrix]
%            .cov        covariance matrices of Gaussian components (must comply with .cov_type)
%                        [dim-by-dim-by-M matrix]
%          X             buffer of N dim-by-1 data set vectors to be evaluated  [dim-by-N]
%          W             (optional) 1-by-N vector of sample weights. If specified, the sample
%                                   set will be weighted according to these weights.
%
%   OUTPUT
%          prior         The prior (without seeing data) probability of a component
%                        density generating any given data vector, i.e. P(C(i)).
%                        This is simply the same as the prior mixing weights,
%                        'gmmDS.weights'. [M-by-1 matrix]
%
%          likelihood    M-by-N martrix where the j,i-th entry is the likelihood
%                        of input column vector i (of X) conditioned on component
%                        density j, i.e. P(X(i)|C(j))
%
%          evidence      1-by-N matrix where the i-th entry is the total data probability
%                        for a given data vector X(i), i.e. P(X(i))=sum_over_all_j[P(X(i)|C(j))]
%
%          posterior     M-by-N matrix where the j,i-th entry is the posterior
%                        probability (after seeing the data) that a component
%                        density j has generated a specific data vector X(i), i.e.
%                        P(C(j)|X(i))   (class posterior probabilities)
%
%
%   See also
%     GMMSAMPLE
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

Nout = nargout;

[dim,nov] = size(X);              % number and dimension of input vectors

if (dim~=gmmDS.dim)
    error(' [ gmmprobability ] Data dimension and GMM model dimension is not the same.');
end

M = gmmDS.M;                      % dumber of component densities
mu    = gmmDS.mu;                 % component means
covar = gmmDS.cov;                % component covariance matrices

prior = gmmDS.weights(:);        % prior mixing probabilities

ones_nov = ones(nov,1);
ones_M   = ones(M,1);

%--- Calculate likelihood
if Nout > 1

  likelihood = zeros(M,nov);        % preallocate component likelihood matrix
  normfact = (2*pi)^(gmmDS.dim/2);  % component density normalizing factor

  switch gmmDS.cov_type             % calculate per component likelihood

  case {'full','diag'}

    for k=1:M,
        cmu = mu(:,k);
        XX = X - cmu(:,ones_nov);
        S = chol(covar(:,:,k))';
        foo = S \ XX;
        likelihood(k,:) = exp(-0.5*sum(foo.*foo, 1))/abs((normfact*prod(diag(S))));
    end

  case {'sqrt','sqrt-diag'}

    for k=1:M,
        cmu = mu(:,k);
        XX = X - cmu(:,ones_nov);
        S = covar(:,:,k);
        foo = S \ XX;
        likelihood(k,:) = exp(-0.5*sum(foo.*foo, 1))/abs((normfact*prod(diag(S))));
    end

  otherwise

    error([' [ gmmprobability ] Unknown covariance type ', mix.cov_type]);

  end

end

likelihood = likelihood + 1e-99;


%--- Calculate evidence
if Nout > 2

  if (nargin == 3)
    evidence = prior' * (likelihood ./ W(ones_M,:));  % weighted
  else
    evidence = prior'*likelihood;                     % non-weighted
  end

  evidence = evidence + 1e-99;

end


%--- Calculate posterior
if Nout > 3

  posterior = likelihood ./ ((1./prior)*evidence) + 1e-99;
  % normalize
  posterior = posterior ./ rvecrep(sum(posterior,1),M);

end



