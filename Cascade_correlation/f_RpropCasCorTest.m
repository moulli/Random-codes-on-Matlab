function [yout, prec] = f_RpropCasCorTest(Xtest, ytest, Wnet, architecture)
% Performs forward prop when given parameters for cascade correlation.
% Xtest is the test set, ytest the labels for the set.
% Wnet is the cell with the matrices provided by the training of net.
% architecture provides in which layer the neurons are.
% yout is the label vector output of the network.
% prec gives the precision of the algorithm based on this yout and ytest.
% Hippolyte MOULLE
    

    %% Initialization:
    Xtest = [ones(1, size(Xtest, 2)); Xtest];
    [n, p] = size(Xtest);
    if size(ytest, 1) > 1
        error("Output vector must be provided as a row vector")
    end
    W0 = Wnet{1};
    
    
    %% Algorithm to build architecture of neural network:
    % Recreating new features:
    Xneur = zeros(n + length(architecture) - 1, p);
    Xneur(1:n, :) = Xtest;
    for i = 2:length(architecture)
        Xneur(n + i - 1, :) = sigmoid(Wnet{i} * Xneur(1:length(Wnet{i}), :));
    end
    % Deducing output:
    output = sigmoid(W0 * Xneur);
    % Recreating label vector:
    [ma, yguess] = max(output);
    yout = yguess - min(yguess) + min(double(ytest));
    
    
    %% Precision of the algorithm:
    prec = sum(yout == ytest) / length(ytest);


end



function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end