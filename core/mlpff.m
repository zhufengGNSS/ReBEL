function Y = mlpff(olType, nodes, X, W1, B1, W2, B2, W3, B3, W4, B4)

% MLPFF  Calculates the output of a ReBEL feed-forward MLP neural network with 2,3 or 4 layers
%
%   Y = mlpff(olType, nodes, X, W1, )
%
%   INPUT
%          olType        output layer type :  'lin' (linear) or 'tanh' (hyperbolic tangent)
%          nodes         neural network layer descriptor vector [num_in num_hid1 (num_hid2) (num_hid3) num_out]
%          X             neural network inputs (xdim-by-M matrix where M is the number of input vectors)
%          Wi            i'th layer weight matrix
%          Bi            i'th layer bias vector
%
%   OUTPUT
%          Y             neural network output (ydim-by-M matrix)
%
%   NOTE : The neural network parameters (weights and biases) can be passed into this function either as seperate
%          unpacked weight and bias matrices (Wi and Bi), or as a single packed vector of parameters (which will be
%          unpacked internally). For the second option 'W1' should be used as the input argument containing the
%          parameter vector.
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

nArgIn = nargin;              % number of input arguments
nLayers = length(nodes)-1;    % number of neural network layers
nInpVecs = size(X,2);         % number of input vectors

if (nArgIn==4)                % unpack weights and biases if needed
    switch nLayers
      case 2
          [W1, B1, W2, B2] = mlpunpack(nodes, W1);
      case 3
          [W1, B1, W2, B2, W3, B3] = mlpunpack(nodes, W1);
      case 4
          [W1, B1, W2, B2, W3, B3, W4, B4] = mlpunpack(nodes, W1);
      otherwise
          error(' [ mlpff ] Only 2, 3 and 4 layer networks supported.');
    end
end


%-- if >= 2 layers

Y = W2*tanh(W1*X + cvecrep(B1,nInpVecs)) + cvecrep(B2,nInpVecs);

%-- if > 2 layers
if (nLayers > 2)
    Y = W3*tanh(Y) + cvecrep(B3,nInpVecs);
end

%-- if > 3 layers
if (nLayers > 3)
    Y = W4*tanh(Y) + cvecrep(B4,nInpVecs);
end


switch olType
  case 'lin'
  case 'tanh'
      Y=tanh(Y);
  otherwise
      error(' [ mlpff ] Output layer type unknown. Only "lin" (linear) and "tanh" (hyperbolic tangent) supported.');
end
