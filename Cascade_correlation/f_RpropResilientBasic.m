function W = f_RpropResilientBasic(X, y, w, lambda, sen)
%% Function that performs resilient backpropagation on a perceptron.
% X is the set, y the labels, w the initial weights matrix.
% theta is the regularization parameter.
% sen is the sensibility of the error that will stop the algorithm.
% Hippolyte MOULLE 


    %%  Initialization:
    W = w;
    % Reconstruction of a usable label vector:
    if size(y, 1) > 1
        error("Output vector must be provided as a row vector")
    end
    class = unique(y);
    yreach = (y == (class(1):class(end))');


    %% Raising errors if dimension or sign problem:
    if size(X, 2) ~= size(yreach, 2)
        error("Not same number of patterns in input and output")
    elseif size(W, 1) ~= size(yreach, 1)
        error("Not same number of weights and classes")
    elseif size(X, 1) ~= size(W, 2)
        error("Different number of features for set and weights")
    elseif lambda < 0
        error("Regularizaton parameter must be positive")
    elseif sen < 0
        error("Sensibility parameter must be positive")
    end

    
    %% Parameters for while loop:
    % Definition of parameters for resilient backprop:
    mu_plus = 1.2;
    mu_moins = 0.5;
    deltamax = 50;
    Delta0 = 0.01 * ones(size(W));
    Loss{1} = zeros(size(W));
    Loss{2} = zeros(size(W));
    % Exit function parameters:
    error1 = 0;
    error2 = 0;
    startCount = 0;
    
    
    %% Main algorithm:
    while error1 - error2 > sen || error1 < error2 || startCount < 2
        
        % Forward propagation:
        yforward = sigmoid(W * X);
        
        % Error computation:
        errorM = yforward - yreach;
        % Scalar error for sensibility criterion:
        errorS = errorM.^2;
        errorS = (1/2) * (1/size(errorM, 2)) * sum(errorS(:));
        error1 = error2;
        error2 = errorS;
        
        % Gradient computation:
        Wreg = [zeros(size(W, 1), 1), W(:, 2:end)]; % regularization
        gradJ = (1/length(y)) * (errorM * X' + lambda * Wreg);
        
        % Resilient backpropagation:
        % Update:
        Loss{1} = Loss{2};
        Loss{2} = gradJ;
        % Sign of derivate multiplication:
        losstemp = Loss{1} .* Loss{2};
        % Delta & which has already reached maximum:
        Delta = mu_plus*Delta0.*(losstemp > 0) + mu_moins*Delta0.*(losstemp < 0) ...
                        + Delta0.*(losstemp == 0);
        Delta(Delta > deltamax) = deltamax;
        Delta0 = Delta;
        % Updating weights:
        DeltaW = -Delta0.*(Loss{2} > 0) + Delta0.*(Loss{2} < 0);
        W = W + DeltaW;
        
        % Loop exit scheme:
        startCount = startCount + 1;
        
    end

    
end



%% Annex sigmoid function:
function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end