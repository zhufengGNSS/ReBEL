% DATAMAT  Packs a vector of data (length N) into a data matrix of dimension M-by-(N-M+1)
%
%   D = datamat(x,M)
%
%   INPUT
%           x       vector of data
%           M       data matrix window (frame) size
%   OUTPUT
%           D       M-by-(N-M+1) datamatrix
%
%
%   Example    D = datamat([1 2 3 4 5 6 7 8 9],3)
%
%     will generate the following datamatrix,
%
%     D = | 3  4  5  6  7  8  9 |
%         | 2  3  4  5  6  7  8 |
%         | 1  2  3  4  5  6  7 |
%
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

%===============================================================================================

function dm=datamat(x,M)

N=length(x);

dm=zeros(N,M);

i=0:N-M;
ii=M:-1:1;
im=repmat(ii',1,N-M+1)+repmat(i,M,1);

dm=x(im);
