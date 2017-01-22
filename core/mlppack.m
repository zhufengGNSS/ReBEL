function [wh, nodes] = mlppack(W1, B1, W2, B2, W3, B3, W4, B4)

% MLPPACK  ReBEL MLP neural network weight matrices vectorizer.
%
%  This function packs the parameters (weights and biases) of ReBEL MLP neural network
%  into a single vector. The network can have between 2, 3 or 4 layers.
%
%   [wh, nodes] = mlppack(W1, B1, W2, B2, W3, B3, W4, B4)
%
% INPUT
%        W1       :   layer 1 weights
%        B1       :   layer 1 biases
%        W2       :   layer 2 weights
%        B2       :   layer 2 biases
%        W3       :   (optional) layer 3 weights
%        B3       :   (optional) layer 3 biases
%        W4       :   (optional) layer 4 weights
%        B4       :   (optional) layer 4 biases
%
% OUTPUT
%        wh       :   vector of 'vectorized' neural network weights
%        nodes    :   MLP neural network layer descriptor vector
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

nIN = nargin;

if (nIN < 4) error(' [ mlppack ] Not enough input arguments. Minimum number of layers is 2.'); end

nLayers = nIN/2;

nodes = zeros(1,nLayers+1);

switch nLayers

 case 2

   [nodes(2) nodes(1)] = size(W1);
   nodes(3) = size(W2,1);
   numParams = nodes(1)*nodes(2) + nodes(2) + nodes(2)*nodes(3) + nodes(3);

   wh=[W1(:) ; B1(:) ; W2(:) ; B2(:)];

 case 3

   [nodes(2) nodes(1)] = size(W1);
   nodes(3) = size(W2,1);
   nodes(4) = size(W3,1);

   numParams = nodes(1)*nodes(2) + nodes(2) + nodes(2)*nodes(3) + nodes(3) + nodes(3)*nodes(4) + nodes(4);

   wh=[W1(:) ; B1(:) ; W2(:) ; B2(:) ; W3(:) ; B3(:)];

 case 4

   [nodes(2) nodes(1)] = size(W1);
   nodes(3) = size(W2,1);
   nodes(4) = size(W3,1);
   nodes(5) = size(W4,1);

   numParams = nodes(1)*nodes(2) + nodes(2) + nodes(2)*nodes(3) + nodes(3) + nodes(3)*nodes(4) + nodes(4) + nodes(4)*nodes(5) + ...
       nodes(5);

   wh=[W1(:) ; B1(:) ; W2(:) ; B2(:) ; W3(:) ; B3(:); W4(:); B4(:)];

 otherwise

  error(' [ mlppack ] MLP neural networks with more than 4 layers are not supported.');

end
