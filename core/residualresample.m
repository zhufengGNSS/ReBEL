function outIndex = residualResample(inIndex, weights);

% RESIDUALRESAMPLE  Residual resampling for SIR. Performs the resampling stage of
%                   the SIR algorithm in order(number of samples) steps.
%
%   outIndex = residualResample(inIndex, weights)
%
%   INPUT
%          inIndex        (r-vector) Input particle indices.
%          weights        (r-vector) Normalised importance weights
%   OUTPUT
%          outIndex          Resampled indices.
%
%   See also
%

%   Copyright  (c) Arnaud Doucet, Nando de Freitas and Rudolph van der Merwe (1998-2002)

%   This file is part of the ReBEL Toolkit. The ReBEL Toolkit is available free for
%   academic use only (see included license file) and can be obtained from
%   http://choosh.csee.ogi.edu/rebel/.  Businesses wishing to obtain a copy of the
%   software should contact rebel@csee.ogi.edu for commercial licensing information.
%
%   See LICENSE (which should be part of the main toolkit distribution) for more
%   detail.

%=============================================================================================

if (nargin ~= 2),
    error(' [ residualResample ] Not enough input arguments.');
end

S = length(weights);       % S = Number of particles.

outIndex = zeros(1,S);   % setup output index buffer

%=== RESIDUAL RESAMPLING  ==========================================================

N_kind= zeros(1,S);

% first integer part
weights_res = S*weights;
N_kind = fix(weights_res);

% residual number of particles to sample
N_res = S-sum(N_kind);

if N_res

    weights_res = (weights_res-N_kind)/N_res;
    cumDist = cumsum(weights_res);

    % generate N_res ordered random variables uniformly distributed in [0,1]
    u = fliplr(cumprod(rand(1,N_res).^(1./(N_res:-1:1))));
    j=1;
    for i=1:N_res
        while (u(1,i)>cumDist(1,j))
            j=j+1;
        end
        N_kind(1,j)=N_kind(1,j)+1;
    end;

end;


%=== COPY RESAMPLED TRAJECTORIES =====================================================

index=1;
for i=1:S
    if (N_kind(1,i)>0)
        for j=index:index+N_kind(1,i)-1
            outIndex(j) = inIndex(i);
        end;
    end;
    index = index+N_kind(1,i);
end


%----------------------------------------------------------------------------------------








