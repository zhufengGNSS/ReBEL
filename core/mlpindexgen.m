function [idxW1, idxB1, idxW2, idxB2, idxW3, idxB3, idxW4, idxB4] = mlpindexgen(nodes)

% MLPINDEXGEN  ReBEL MLP neural network parameter matrices de-vectorizing index generator
%
%  This function generates the needed index vectors to directly devoctorize a single ReBEL MLP
%  neural network parameter vector into the corresponding weight and bias matrices. The output
%  arguments are the index vectors for each layers weight and bias matrices. 'nodes' specify
%  the MLP structure.
%
%  [idxW1, idxB1, idxW2, idxB2, idxW3, idxB3, idxW4, idxB4] = mlpindexgen(nodes)
%
%  INPUT
%        nodes    :   MLP neural network layer descriptor vector
%
%  OUTPUT
%        idxW1       :   layer 1 weights index vector
%        idxB1       :   layer 1 biases index vector
%        idxW2       :   layer 2 weights index vector
%        idxB2       :   layer 2 biases index vector
%        idxW3       :   (optional) layer 3 weights index vector
%        idxB3       :   (optional) layer 3 biases index vector
%        idxW4       :   (optional) layer 4 weights index vector
%        idxB4       :   (optional) layer 4 biases index vector
%
%
%  SEE ALSO:
%            mlppack, mlpunpack
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


nLayers = length(nodes)-1;

if (nLayers<2)
  error(' [ mlpindexgen ]  MLP neural networks need at least 2 layers.');
elseif (nLayers>4)
  error(' [ mlpindexgen ]  MLP neural networks with more than 4 layers are not supported.');
end


%--- If nLayers at leat == 2

  numW1 = nodes(1)*nodes(2);        % number of parameters in W1 matrix
  numB1 = nodes(2);                 % number of parameters in B1 matrix (actually a vector)
  numW2 = nodes(2)*nodes(3);        % number of parameters in W2
  numB2 = nodes(3);                 % number of parameters in B2

  i=0;
  j=i+numW1; idxW1 = i+1:j; i=j;
  j=i+numB1; idxB1 = i+1:j; i=j;
  j=i+numW2; idxW2 = i+1:j; i=j;
  j=i+numB2; idxB2 = i+1:j; i=j;


%--- If nLayers at least == 3
if (nLayers > 2)

   numW3 = nodes(3)*nodes(4);
   numB3 = nodes(4);

   j=i+numW3; idxW3 = i+1:j; i=j;
   j=i+numB3; idxB3 = i+1:j; i=j;

end


%--- If nLayers == 4
if (nLayers > 3)

   numW4 = nodes(4)*nodes(5);
   numB4 = nodes(5);

   j=i+numW4; idxW4 = i+1:j; i=j;
   j=i+numB4; idxB4 = i+1:j; i=j;

end
