function W = f_PolakResilientBasic(X, y, w, lambda, sen)
%% Function that performs Polak-Ribière optimisation on a perceptron.
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
    % Definition of parameters for Polak-Ribière:
    errorM = sigmoid(W * X) - yreach;
    Wreg = [zeros(size(W, 1), 1), W(:, 2:end)]; % regularization
    Grad{2} = (1/length(y)) * (errorM * X' + lambda * Wreg);
    Grad{1} = zeros(size(W));
    dk = Grad{2};
    startCount = 0;
    
    
    %% Main algorithm:
    while sqrt(sum(sum(Grad{2}.^2))) >= sen
        
        % INITIALISER AVEC UN GRADIENT SIMPLE
        
        % Armijo for the step:
        alphak = 1;
        tau = 0.01;
        omega1 = 0.0001;
        omega2 = 0.99;
        Wreg = [zeros(size(W, 1), 1), W(:, 2:end)];
        errorM2 = sigmoid((W+alphak*dk) * X) - yreach;
        GradM2 = (1/length(y)) * (errorM2 * X' + lambda * Wreg);
        sum(sum(Grad{2}.*dk))
        sum(sum(errorM.^2))
        while sum(sum(errorM2.^2)) > sum(sum(errorM.^2)) + omega1*alphak*sum(sum(Grad{2}.*dk)) ...
                || sum(sum(GradM2.*dk)) < omega2*sum(sum(Grad{2}.*dk))
            alphak = alphak * ((1 - 2*tau) * (rand - 0.5) + 0.5);
            % Wolfe criterion:
            errorM2 = sigmoid((W+alphak*dk) * X) - yreach;
            GradM2 = (1/length(y)) * (errorM2 * X' + lambda * Wreg);
        end
        
        % Weight and gradient update:
        W = W + alphak*dk;
        Grad{1} = Grad{2};
        errorM = sigmoid(W * X) - yreach;
        Wreg = [zeros(size(W, 1), 1), W(:, 2:end)];
        Grad{2} = (1/length(y)) * (errorM * X' + lambda * Wreg);
        
        % Descent direction update:
        bk = sum(sum(Grad{2}.*(Grad{2} - Grad{1}))) / sum(sum(Grad{1}.*Grad{1}));
        dk = -Grad{2} + bk * dk;
        
        % Loop exit scheme:
        startCount = startCount + 1
        
    end

    
end



%% Annex sigmoid function:
function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end