function NoiseDS = convgausns(NoiseDS, target_cov_type)

% CONVGAUSNS  Convert Gaussian noise source from one cov_type to another
%
%   NoiseDS = convgausns(NoiseDS, target_cov_type)
%
%   Input
%          NoiseDS        : (NoiseDS) Input noise source data structure (this must be of type 'gaussian' or 'combo-gaussian')
%          target_cov_type : (string)  Target cov_type : 'full', 'diag', 'sqrt', 'sqrt-diag'
%
%   Output
%          NoiseDS        : (NoiseDS) Converted noise source data structure
%
%
%   See also
%     GENNOISEDS
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


%--- ERROR CHECKING ---------------------------------------------------------------------

if (nargin < 2) error(' [ convgausns ] Incorrect number of input arguments.'); end

error(consistent(NoiseDS,'NoiseDS'));       % Check general consistency of NoiseDS data structure

if ~stringmatch(NoiseDS.ns_type,{'gaussian','combo-gaussian','gmm'})
  error(' [ convgausns ] This function can only operate on Gaussian or Combination-Gaussian noise sources.');
end

if ~ischar(target_cov_type)
    error(' [ convgausns ] Second input argument must be a string.');
end

if ~stringmatch(target_cov_type,{'full','diag','sqrt','sqrt-diag'})
    error([' [ convgausns ] Unknown target cov_type ''' target_cov_type '''']);
end


%----------------------------------------------------------------------------------------
switch target_cov_type

%........................................................................................
case 'full'

  switch NoiseDS.cov_type
  case {'sqrt','sqrt-diag'}
      if stringmatch(NoiseDS.ns_type,'gmm')
          for k=1:NoiseDS.M,
              NoiseDS.cov(:,:,k) = NoiseDS.cov(:,:,k) * NoiseDS.cov(:,:,k)';
          end
      else
          NoiseDS.cov = NoiseDS.cov * NoiseDS.cov';
      end
  end
  NoiseDS.cov_type = target_cov_type;

%........................................................................................
case 'diag'

  switch NoiseDS.cov_type
  case {'sqrt','sqrt-diag'}
      if stringmatch(NoiseDS.ns_type,'gmm')
          for k=1:NoiseDS.M,
              NoiseDS.cov(:,:,k) = diag(diag(NoiseDS.cov(:,:,k)*NoiseDS.cov(:,:,k)'));
          end
      else
          NoiseDS.cov = diag(diag(NoiseDS.cov*NoiseDS.cov'));
      end
  otherwise
      if stringmatch(NoiseDS.ns_type,'gmm')
          for k=1:NoiseDS.M,
              NoiseDS.cov(:,:,k) = diag(diag(NoiseDS.cov(:,:,k)));
          end
      else
          NoiseDS.cov = diag(diag(NoiseDS.cov));
      end
  end
  NoiseDS.cov_type = target_cov_type;

%........................................................................................
case 'sqrt'

  switch NoiseDS.cov_type
  case {'full','diag'}
      if stringmatch(NoiseDS.ns_type,'gmm')
          for k=1:NoiseDS.M,
              NoiseDS.cov(:,:,k) = chol(NoiseDS.cov(:,:,k))';
          end
      else
          NoiseDS.cov = chol(NoiseDS.cov)';
      end
  end
  NoiseDS.cov_type = target_cov_type;

%........................................................................................
case 'sqrt-diag'

  switch NoiseDS.cov_type
  case {'full','diag'}
      if stringmatch(NoiseDS.ns_type,'gmm')
          for k=1:NoiseDS.M,
              NoiseDS.cov(:,:,k) = diag(diag(chol(NoiseDS.cov(:,:,k))'));
          end
      else
          NoiseDS.cov = diag(diag(chol(NoiseDS.cov)'));
      end
   otherwise
      if stringmatch(NoiseDS.ns_type,'gmm')
          for k=1:NoiseDS.M,
              NoiseDS.cov(:,:,k) = diag(diag(NoiseDS.cov(:,:,k)));
          end
      else
          NoiseDS.cov = diag(diag(NoiseDS.cov));
      end
  end
  NoiseDS.cov_type = target_cov_type;

%........................................................................................
otherwise
  error([' [ convgausns ] Unknown target cov_type ''' target_cov_type '''']);
end


