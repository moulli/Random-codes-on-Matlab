function out = fReinforcement(Winit, Xcell, step, lambda)
% Function that performs resilient back propagation for reinforcement
% learning on a neural net.
% Winit is the weight we want to update
% Xdata must be provided as {s x 2} cell, with:
% - s the number of trials we ran (each made of several actions)
%   for each trial, data are line vectors one after the other
% - one column for the data, and another one +1 or -1 depending on result
% reward is higher is we want to highly reward the program 
% lambda is the regularization parameter
% out is a cell containing the updated weights


    % Creating usable dataset:
    [setn, ~] = size(Xcell);
    Xdata = [];
    yreward = [];
    for i = 1:setn
        Xdata = [Xdata, Xcell{i, 1}'];
        yreward = [yreward, Xcell{i, 2} * ones(1, size(Xcell{i, 1}, 1))];
    end
    
    % Info & adding bias:
    [~, xm] = size(Xdata);
    Xdata = [ones(1, xm); Xdata];
    
    % Neural net structure:
    nlayer = length(Winit) - 1;
    W = cell(nlayer + 1, 2);
    W{1, 1} = Winit{1}; % Size of weight matrices
    W{1, 2} = [1, numel(W{1, 1})]; % First and last element in the unroll
    for i = 2:nlayer
        W{i, 1} = Winit{i};
        W{i, 2} = [W{i-1, 2}(2)+1, W{i-1, 2}(2)+numel(W{i, 1})];
    end
    W{nlayer+1, 1} = Winit{nlayer+1};
    W{nlayer+1, 2} = [W{nlayer, 2}(2)+1, W{nlayer, 2}(2)+numel(W{nlayer+1, 1})];
    % Creating a big unrolled matrix to manipulate the data:
    Wlayed = zeros(W{nlayer+1, 2}(2), 1);
    for i = 1:nlayer+1
        Wlayed(W{i, 2}(1):W{i, 2}(2)) = reshape(W{i, 1}, [numel(W{i, 1}), 1]);
    end

    % Keeping track of each layer result:
    A = cell(nlayer+1, 1);
    error = cell(nlayer+1, 1);

    % Forward propagation:
    Z1 = W{1, 1} * Xdata;
    A1 = sigmoid(Z1);
    A{1} = [ones(1, size(A1, 2)); A1];
    for j = 2:nlayer
        Zj = W{j, 1} * A{j-1};
        Aj = sigmoid(Zj);
        A{j} = [ones(1, size(Aj, 2)); Aj];
    end
    Znlayer1 = W{nlayer+1, 1} * A{nlayer};
    A{nlayer+1} = sigmoid(Znlayer1);
    yforward = A{nlayer+1};

    % Computation of error for reinforcement learning:
    [~, ykeep] = max(yforward);
    yerror = yreward .* log(ykeep .* yforward);

    % Computation of derivatives:
    error = cell(nlayer+1, 1);
    gradJ = cell(nlayer+1, 1);
    error{nlayer+1} = yerror;
    Wreg = [zeros(size(W{nlayer+1, 1}, 1), 1), W{nlayer+1, 1}(:, 2:end)];
    gradJ{nlayer+1} = (1 / xm) * (error{nlayer+1} * A{nlayer}' + lambda * Wreg); 
    for j = fliplr(2:nlayer)
        errorj = (W{j+1, 1}' * error{j+1}) .* A{j} .* (1 - A{j});
        error{j} = errorj(2:end, :);
        Wreg = [zeros(size(W{j, 1}, 1), 1), W{j, 1}(:, 2:end)];
        gradJ{j} = (1 / xm) * (error{j} * A{j-1}' + lambda * Wreg); 
    end
    error1 = (W{2, 1}' * error{2}) .* A{1} .* (1 - A{1});
    error{1} = error1(2:end, :);
    Wreg = [zeros(size(W{1, 1}, 1), 1), W{1, 1}(:, 2:end)];
    gradJ{1} = (1 / xm) * (error{1} * Xdata' + lambda * Wreg);

    % Total gradient computation:
    gradT = zeros(size(Wlayed));
    for j = 1:(nlayer+1)
        gradT(W{j, 2}(1):W{j, 2}(2)) = gradJ{j}(:);
    end

    % Update:
    Wlayed = Wlayed - step * gradT;
    
    % We return the weights:
    out = cell(nlayer+1, 1);
    for i = 1:(nlayer+1)
        out{i} = reshape(Wlayed(W{i, 2}(1):W{i, 2}(2)), size(W{i, 1}));
    end

end



function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end