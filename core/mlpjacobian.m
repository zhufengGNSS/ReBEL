function [J1,J2] = mlpjacobian(jacType, olType, nodes, X, W1, B1, W2, B2, W3, B3, W4, B4)

% MLPJACOBIAN   Calculates the Jacobian (first partial derivative matrix)
%               of a ReBEL MLP neural network. The independent variable can
%               be either the network input 'X' or the network parameters (weights
%               and biases packed into a single vector).
%
%  [J1,J2] = mlpjacobian(jacType, olType, nodes, X, W1, B1, W2, B2, W3, B3, W4, B4)
%
%  INPUT
%         jacType :   Jacobian type : 'dydx'  : d(output)/d(input)
%                                      'dydw'  : d(output)/d(parameters)
%                                      'dydxw' : J1=d(output)/d(input) J2=d(output)/d(parameters)
%         olType  :   output layer type, linear ('lin') or hyperbolic-tangent ('tanh')
%         nodes    :   network layer descriptor vector  [numIn numHid1 (numHid2) (numHid3) numOut]
%         X        :   neural network input
%         W1       :   layer 1 weights
%         B1       :   layer 1 biases
%         W2       :   layer 2 weights
%         B2       :   layer 2 biases
%         W3       :   (optional) layer 3 weights
%         B3       :   (optional) layer 3 biases
%         W4       :   (optional) layer 4 weights
%         B4       :   (optional) layer 4 biases
%
%         NOTE  - If only W1 is specified, then it is assumed that W1 is a packed vector containing ALL
%                 the neural network parameters (weights and biases for all layers).
%
%  OUTPUT
%         J1       :   Jacobian matrix: This matrix has dimensions :
%                        (dimension of network output)-by-(dimension of independent variable)
%         J2       :   (optional) if jacType='dydxw' then J1=d(out)/d(X) and J2=d(out)/d(parameters)
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

%==================================================================================================

numArgIn = nargin;

if (numArgIn == 5),

  wh = W1;
  numWeights = length(wh);
  numLayers = length(nodes)-1;
  switch numLayers
   case 2
     [W1,B1,W2,B2] = mlpunpack(nodes, wh);
   case 3
     [W1,B1,W2,B2,W3,B3] = mlpunpack(nodes, wh);
   case 4
     [W1,B1,W2,B2,W3,B3,W4,B4] = mlpunpack(nodes, wh);
   otherwise
     error(' [ mlpjacobian ] ReBEL MLP neural networks can only have 2, 3 or 4 weight layers.');
   end

else

   numWeights = sum([nodes(1:end-1).*nodes(2:end) nodes(2:end)]);
   numLayers  = (numArgIn - 4)/2;

end


hiddenLayer1Act = tanh(W1*X + B1);                     % first hidden layer output;


switch numLayers
case 2
  outputLayerAct = W2*hiddenLayer1Act + B2;            % calculate output layer activation
  nOUT = length(B2);                                   % number of output units
case 3
  hiddenLayer2Act = tanh(W2*hiddenLayer1Act + B2);     % calculate second hidden layer output
  outputLayerAct  = W3*hiddenLayer2Act + B3;           % calculate output layer activation
  nOUT = length(B3);                                   % number of output units
case 4
  hiddenLayer2Act = tanh(W2*hiddenLayer1Act + B2);     % calculate second hidden layer output
  hiddenLayer3Act = tanh(W3*hiddenLayer2Act + B3);     % calculate third hidden layer output
  outputLayerAct  = W4*hiddenLayer3Act + B4;           % calculate output layer activation
  nOUT = length(B4);                                   % number of output units
end

% Output layer
switch olType
case 'lin'
case 'tanh'
  outputLayerAct  = tanh(outputLayerAct);               % calculate output layer output
otherwise
  error(' [ mlpjacobian ] Unknown output layer activation function');
end


% Deltas for output layer
switch (olType),
 case 'lin'
  deltaOUT = eye(nOUT);
 case 'tanh'
  deltaOUT = diag(1-outputLayerAct.^2);
end

% Deltas for hidden layers
if (nOUT>1)
    switch numLayers
    case 2
        delta1 = cvecrep(1-hiddenLayer1Act.^2, nOUT) .* (W2'*deltaOUT);
    case 3
        delta2 = cvecrep(1-hiddenLayer2Act.^2, nOUT) .* (W3'*deltaOUT);
        delta1 = cvecrep(1-hiddenLayer1Act.^2, nOUT) .* (W2'*delta2);
    case 4
        delta3 = cvecrep(1-hiddenLayer3Act.^2, nOUT) .* (W4'*deltaOUT);
        delta2 = cvecrep(1-hiddenLayer2Act.^2, nOUT) .* (W3'*delta3);
        delta1 = cvecrep(1-hiddenLayer1Act.^2, nOUT) .* (W2'*delta2);
    end
else
    switch numLayers
    case 2
        delta1 = (1-hiddenLayer1Act.^2) .* (W2'*deltaOUT);
    case 3
        delta2 = (1-hiddenLayer2Act.^2) .* (W3'*deltaOUT);
        delta1 = (1-hiddenLayer1Act.^2) .* (W2'*delta2);
    case 4
        delta3 = (1-hiddenLayer3Act.^2) .* (W4'*deltaOUT);
        delta2 = (1-hiddenLayer2Act.^2) .* (W3'*delta3);
        delta1 = (1-hiddenLayer1Act.^2) .* (W2'*delta2);
    end
end



% Calculate the apropriate Jacobian

switch (jacType)

%--- Derivative with respect to input
case 'dydx'

  J1 = (W1' * delta1)';

%--- Derivative with respect to parameters
case {'dydw','dydxw'},

  dW = zeros(numWeights,nOUT);

  for j=1:nOUT,

    switch numLayers
    case 2
      ddB1 = delta1(:,j);
      ddW1 = ddB1*X';
      ddB2 = deltaOUT(:,j);
      ddW2 = ddB2*hiddenLayer1Act';
      dW(:,j) = [ddW1(:) ; ddB1 ; ddW2(:) ; ddB2];
    case 3
      ddB1 = delta1(:,j);
      ddW1 = ddB1*X';
      ddB2 = delta2(:,j);
      ddW2 = ddB2*hiddenLayer1Act';
      ddB3 = deltaOUT(:,j);
      ddW3 = ddB3*hiddenLayer2Act';
      dW(:,j) = [ddW1(:) ; ddB1 ; ddW2(:) ; ddB2 ; ddW3(:) ; ddB3];
    case 4
      ddB1 = delta1(:,j);
      ddW1 = ddB1*X';
      ddB2 = delta2(:,j);
      ddW2 = ddB2*hiddenLayer1Act';
      ddB3 = delta3(:,j);
      ddW3 = ddB3*hiddenLayer2Act';
      ddB4 = deltaOUT(:,j);
      ddW4 = ddB4*hiddenLayer3Act';
      dW(:,j) = [ddW1(:) ; ddB1 ; ddW2(:) ; ddB2 ; ddW3(:) ; ddB3 ; ddW4(:) ; ddB4];
    end

  end

  switch jacType
  case  'dydxw'
    J1 = (W1' * delta1)';
    J2 = dW';
  otherwise
    J1 = dW';
  end

  % test by perturbation
  %if (strcmp(deriv_type,'dydxp')),
  %  epsilon=1e-8;
  %  WW = cvecrep(W,length(W)) + epsilon*eye(length(W));
  %  YY = zeros(length(W),1);
  %  dYdW=zeros(length(W),1);
  %  Y  = lo{end};
  %  for j=1:length(W),
  %    YY(j) = nnetN(olType,WW(:,j),nodes,X);
  %    dYdW(j) = (YY(j)-Y)/epsilon;
  %  end
  %  ZZ=(J1'-dYdW)./(dYdW);
  %  disp(['Max discrepancy : ' num2str(max(abs(ZZ)))]);
  %end

otherwise

  error(' [ mlpjacobian ] Unknown Jacobian type.');

end


if (0)


  % test by perturbation
  switch (jacType)

  case 'dydx'
    epsilon=1e-8;
    XX = cvecrep(X,length(X)) + epsilon*eye(length(X));
    YY = zeros(nOUT,length(X));
    dYdX=zeros(nOUT,length(X));
    Y  = outputLayerAct;
    for j=1:length(X),
      YY(:,j) = nnetN(olType,wh,nodes,XX(:,j));
      dYdW(:,j) = (YY(:,j)-Y)/epsilon;
    end
    ZZ=(J1-dYdW)./(dYdW);
    disp(['Max discrepancy : ' num2str(max(max(abs(ZZ))))]);

  case 'dydw'
    epsilon=1e-8;
    WW = cvecrep(wh,length(wh)) + epsilon*eye(length(wh));
    YY = zeros(nOUT,length(wh));
    dYdW=zeros(nOUT,length(wh));
    Y  = outputLayerAct;
    for j=1:length(wh),
      YY(:,j) = nnetN(olType,WW(:,j),nodes,X);
      dYdW(:,j) = (YY(:,j)-Y)/epsilon;
    end
    ZZ=(J1-dYdW)./(dYdW);
    disp(['Max discrepancy : ' num2str(max(max(abs(ZZ))))]);

  end

end
