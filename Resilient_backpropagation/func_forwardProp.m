function out = func_forwardProp(W, Xtest)
% Function that will simply perform forward propagation.
% Returns vector with associated class for each example.

    % Getting info:
    tries = length(W);
    
    % Transforming data to use them:
    Xtest = [ones(1, size(Xtest, 2)); Xtest];
    
    % Forward propagation:
    Aj = sigmoid(W{1}*Xtest);
    Aj = [ones(1, size(Aj, 2)); Aj];
    for j = 2:(tries-1)
        Aj = sigmoid(W{j}*Aj);
        Aj = [ones(1, size(Aj, 2)); Aj];
    end
    Aj = sigmoid(W{tries}*Aj);
    yforward = Aj;
    
    % Transforming into usable data:
    [ymax, out] = max(yforward);

end



function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end