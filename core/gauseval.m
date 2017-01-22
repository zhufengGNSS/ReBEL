function likelihood = gauseval(gausDS, X, logflag)

% GAUSEVAL  Calculates the likelihood of a dataset X under a Gaussian probability
%           density described by the Gaussian data structure 'gausDS',
%           i.e. P(X|gausDS). The column vectors of X are treated as IID samples to be
%           evaluated. The function return a likelihood row-vector that contains one
%           likelihood value for each IID sample in X.
%
%   likelihood = gauseval(gausDS, X)
%
%   INPUT
%          gausDS       Gaussian data structure with the following fields
%             .cov_type (string) covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'
%             .dim      (scalar)   dimension
%             .mu       (c-vector) mean vector  (dim-by-1)
%             .cov      (matrix)   covariance matrix of type cov_type  (dim-by-dim)
%          X            (dim-by-M) buffer of M dim-by-1 data set vectors to be evaluated
%          logflag      <optional> if 'logflag==1' then the probability is calculated in the log
%                                  domain, i.e. likelihood=log(likelihood) (default : logflag=0)
%   OUTPUT
%          likelihood   (1-by-M) buffer of likelihood values
%
%   See also
%     GAUSSAMP
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

switch nargin
case 3
case 2
  logflag=0;
otherwise
  error(' [ gauseval ] Not enough inputs.');
end


nov = size(X, 2);                % number of input vectors
likelihood = zeros(1, nov);      % preallocate likelihood matrix


switch gausDS.cov_type           % calculations depend on covariance type

case {'full','diag'}

    normfact = (2*pi)^(gausDS.dim/2);
    XX = X - cvecrep(gausDS.mu, nov);
    S = chol(gausDS.cov)';
    foo = S \ XX;
    if logflag
        likelihood = -0.5*sum(foo.*foo, 1) - log(normfact*abs(prod(diag(S))));
    else
        likelihood = exp(-0.5*sum(foo.*foo, 1))./(normfact*abs(prod(diag(S))));
    end

case {'sqrt','sqrt-diag'}

    normfact = (2*pi)^(gausDS.dim/2);
    XX = X - cvecrep(gausDS.mu, nov);
    S = gausDS.cov;
    foo = S \ XX;
    if logflag
        likelihood = -0.5*sum(foo.*foo, 1) - log(normfact*abs(prod(diag(S))));
    else
        likelihood = exp(-0.5*sum(foo.*foo, 1))./(normfact*abs(prod(diag(S))));
    end


otherwise

    error([' [ gauseval ] Unknown covariance type ', mix.cov_type]);

end

