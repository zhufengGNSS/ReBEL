function [W1, B1, W2, B2, W3, B3, W4, B4] = mlpweightinit(nodes)

% MLPWEIGHTINIT   Initializes the weights of a ReBEL MLP feedforward neural network
%
%    [W1, B1, W2, B2, W3, B3, W4, B4] = mlpweightinit(nodes)
%
%    INPUT
%             nodes        : (r-vector)  neural network layer descriptor vector [num_in num_hid1 (num_hid2) (num_hid3) num_out]
%
%    OUTPUT
%             Wi           : (matrix) ith layer weight matrix
%             Bi           : (matrix) ith layer bias vector
%
%       EXAMPLE: [W1,B1,W2,B2] = mlpweightinit([3 5 2])
%            Returns the weights for a standard 2 layer network
%            with 3 inputs, 5 hidden units, and 2 outputs.
%
%       EXAMPLE: [W1,B1,W2,B2,W3,B3] = mlpweightinit([4 5 5 2])
%            Returns the weights of a 3 layer network with 4 inputs,
%            5 units in the first hidden layer, 5 units in the second
%            hidden layer and 2 output units.
%
%       EXAMPLE: [W1,B1,W2,B2,W3,B3,W4,B4] = mlpweightinit([5 4 2 4 1])
%            Returns the weights of a 4 layer network with 5 inputs,
%            4 units in the first hidden layer, 2 units in the second
%            hidden layer, 5 units in the third hidden layer and 1 output unit.
%
%   NOTE : If only one output argument is given, i.e. W = mlpWeightInit(nodes),
%          then all the neural network paramaters (weights & biases) are returned
%          in a packed (single vector) format. See 'mlppack' and 'mlpunpack'
%
%
%   See also
%   MLPPACK, MLPUNPACK, MLPINDEXGEN
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

rand('state',sum(100*clock));   % stir the pot a bit...

lenNodes = length(nodes);


W1 = (rand(nodes(2),nodes(1))-0.5) *sqrt(3/nodes(1));
B1 = (rand(nodes(2),1)-0.5) *sqrt(3/nodes(1));


% two layers
if lenNodes > 2,
    W2 = (rand(nodes(3),nodes(2))-0.5) *sqrt(3/nodes(2));
    B2 = (rand(nodes(3),1)-0.5) *sqrt(3/nodes(2));
end

% three layers
if lenNodes > 3,
    W3 = (rand(nodes(4),nodes(3))-0.5) *sqrt(3/nodes(3));
    B3 = (rand(nodes(4),1)-0.5) *sqrt(3/nodes(3));
end

% four layers
if lenNodes > 4
    W4 = (rand(nodes(5),nodes(4))-0.5) *sqrt(3/nodes(4));
    B4 = (rand(nodes(5),1)-0.5) *sqrt(3/nodes(4));
end

if lenNodes > 5
    error(' [ mlpWeightInit ] Only 2, 3 and 4 layer MLp neural nets supported.');
end


%--- If requesting packed weights for output

if (nargout < 2),

  switch lenNodes

   case 3
     W1 = mlppack(W1,B1,W2,B2);

   case 4
     W1 = mlppack(W1,B1,W2,B2,W3,B3);

   case 5
     W1 = mlppack(W1,B1,W2,B2,W3,B3,W4,B4);

  end

end

