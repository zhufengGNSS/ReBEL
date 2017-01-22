function dups = checkdups(x)

% CHECKDUPS  Checks for the presence of duplicate entries in an array
%
%     dups = CHECKDUPS(x)
%
%     dups = number_of_duplicates if there are duplicate entries in x
%     otherwise dups = 0

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

x=x(:);

if (length(x)>1)

  x = sort(x);
  y = [x(end); x(1:end-1)];
  d = x-y;

  dups = sum(d==0);

else

  dups = 0;

end

