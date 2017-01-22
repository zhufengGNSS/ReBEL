function X = gaussample(gausDS, N)

% GAUSSAMPLE  Draw N samples from the Gaussian distribution (pdf) described by the
%             Gaussian data structure 'gausDS'.
%
%   X = gaussample(gaussDS, X)
%
%   INPUT
%          gausDS       Gaussian data structure with the following fields
%             .cov_type (string)   covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'
%             .dim      (scalar)   dimension
%             .mu       (c-vector) mean vector  (dim-by-1)
%             .cov      (matrix)   covariance matrix of type cov_type  (dim-by-dim)
%          N            (scalar)   number of samples to generate
%   OUTPUT
%          X            (matrix)   buffer of generated samples (dim-by-N)
%
%   See also
%     GAUSEVAL
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

switch gausDS.cov_type           % calculations depend on covariance type
case {'full','diag'}
    S = chol(gausDS.cov)';
case {'sqrt','sqrt-diag'}
    S = gausDS.cov;
otherwise
    error([' [ gaussample ] Unknown covariance type ', mix.cov_type]);
end

X = S * randn(gausDS.dim,N) + gausDS.mu(:,ones(N,1));