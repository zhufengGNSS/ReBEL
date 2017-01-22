function [W1, B1, W2, B2, W3, B3, W4, B4] = mlpunpack(nodes, wh)

% MLPUNPACK  ReBEL MLP neural network weight matrices de-vectorizer.
%
%  This function unpacks the parameters (weights and biases) of ReBEL MLP neural network
%  from a single vector into the correct weight and bias matrices as specified by the
%  neural network layer descriptor vector, nodes. Only 2,3 and 4 layer networks are supported.
%
%   [W1, B1, W2, B2, W3, B3, W4, B4] = mlpunpack(nodes, wh)
%
%  INPUT
%        wh       :   vector of 'vectorized' neural network weights (created with 'nnpack.m')
%        nodes    :   MLP neural network layer descriptor vector
%
%  OUTPUT
%        W1       :   layer 1 weights
%        B1       :   layer 1 biases
%        W2       :   layer 2 weights
%        B2       :   layer 2 biases
%        W3       :   (optional) layer 3 weights
%        B3       :   (optional) layer 3 biases
%        W4       :   (optional) layer 4 weights
%        B4       :   (optional) layer 4 biases
%
%
%  SEE ALSO:
%            mlppack
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

if (nargin ~= 2) error(' [ mlpunpack ] Not enough input arguments.'); end

nLayers = length(nodes)-1;


switch nLayers

 case 2

   numW1 = nodes(1)*nodes(2);
   numB1 = nodes(2);
   numW2 = nodes(2)*nodes(3);
   numB2 = nodes(3);

   W1=zeros(nodes(2),nodes(1));
   B1=zeros(nodes(2),1);
   W2=zeros(nodes(3),nodes(2));
   B2=zeros(nodes(3),1);

   i=0;
   j=i+numW1; W1(:) = wh(i+1:j); i=j;
   j=i+numB1; B1 = wh(i+1:j); i=j;
   j=i+numW2; W2(:) = wh(i+1:j); i=j;
   j=i+numB2; B2 = wh(i+1:j);

 case 3

   numW1 = nodes(1)*nodes(2);
   numB1 = nodes(2);
   numW2 = nodes(2)*nodes(3);
   numB2 = nodes(3);
   numW3 = nodes(3)*nodes(4);
   numB3 = nodes(4);

   W1=zeros(nodes(2),nodes(1));
   B1=zeros(nodes(2),1);
   W2=zeros(nodes(3),nodes(2));
   B2=zeros(nodes(3),1);
   W3=zeros(nodes(4),nodes(3));
   B3=zeros(nodes(4),1);

   i=0;
   j=i+numW1; W1(:) = wh(i+1:j); i=j;
   j=i+numB1; B1 = wh(i+1:j); i=j;
   j=i+numW2; W2(:) = wh(i+1:j); i=j;
   j=i+numB2; B2 = wh(i+1:j); i=j;
   j=i+numW3; W3(:) = wh(i+1:j); i=j;
   j=i+numB3; B3 = wh(i+1:j);

 case 4

   numW1 = nodes(1)*nodes(2);
   numB1 = nodes(2);
   numW2 = nodes(2)*nodes(3);
   numB2 = nodes(3);
   numW3 = nodes(3)*nodes(4);
   numB3 = nodes(4);
   numW4 = nodes(4)*nodes(5);
   numB4 = nodes(5);

   W1=zeros(nodes(2),nodes(1));
   B1=zeros(nodes(2),1);
   W2=zeros(nodes(3),nodes(2));
   B2=zeros(nodes(3),1);
   W3=zeros(nodes(4),nodes(3));
   B3=zeros(nodes(4),1);
   W4=zeros(nodes(5),nodes(4));
   B4=zeros(nodes(5),1);

   i=0;
   j=i+numW1; W1(:) = wh(i+1:j); i=j;
   j=i+numB1; B1 = wh(i+1:j); i=j;
   j=i+numW2; W2(:) = wh(i+1:j); i=j;
   j=i+numB2; B2 = wh(i+1:j); i=j;
   j=i+numW3; W3(:) = wh(i+1:j); i=j;
   j=i+numB3; B3 = wh(i+1:j);
   j=i+numW4; W4(:) = wh(i+1:j); i=j;
   j=i+numB4; B4 = wh(i+1:j);


 otherwise

  error(' [ mlpunpack ] MLP neural networks with more than 4 layers are not supported.');

end
