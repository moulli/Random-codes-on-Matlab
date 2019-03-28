function out = func_neuralNet(Xtrain, ytrain, structure, lambda, epoch)
% Function that performs resilient back propagation on a neural net.

    % Info & adding bias:
    [xn, xm] = size(Xtrain);
    Xtrain = [ones(1, size(Xtrain, 2)); Xtrain];
    nmov = length(unique(ytrain));
    ytemp = zeros(nmov, length(ytrain));
    for j = 1:length(ytrain)
        for i = 1:nmov
            ytemp(i, j) = (i == ytrain(j));
        end
    end
    ytrain = ytemp;
    
    % Neural net structure:
    nlayer = length(structure);
    W = cell(nlayer+1, 2);
    W{1, 1} = [structure(1), xn+1]; % Size of weight matrices
    W{1, 2} = [1, prod(W{1, 1})]; % First and last element in the unroll
    for i = 2:nlayer
        W{i, 1} = [structure(i), structure(i-1)+1];
        W{i, 2} = [W{i-1, 2}(2)+1, W{i-1, 2}(2)+prod(W{i, 1})];
    end
    W{nlayer+1, 1} = [nmov, structure(nlayer)+1];
    W{nlayer+1, 2} = [W{nlayer, 2}(2)+1, W{nlayer, 2}(2)+prod(W{nlayer+1, 1})];
    % Creating a big unrolled matrix to manipulate the data:
    Wlayed = 0.012 * (rand(W{nlayer+1, 2}(2), 1) - 0.5);
    
    % Definition of parameters for resilient backprop:
    mu_plus = 1.2;
    mu_moins = 0.5;
    deltamax = 50;
    Delta0 = 0.01 * ones(size(Wlayed));
    Loss{1} = zeros(size(Wlayed));
    Loss{2} = zeros(size(Wlayed));
    
    % Optimization using resilient backpropagation:
    for i = 1:epoch

        % Keeping track of each layer result:
        A = cell(nlayer+1, 1);
        Weight = cell(nlayer+1, 1);
        
        % Forward propagation:
        Weight{1} = reshape(Wlayed(W{1, 2}(1):W{1, 2}(2)), W{1, 1});
        Z1 = Weight{1} * Xtrain;
        A1 = sigmoid(Z1);
        A{1} = [ones(1, size(A1, 2)); A1];
        for j = 2:nlayer
            Weight{j} = reshape(Wlayed(W{j, 2}(1):W{j, 2}(2)), W{j, 1});
            Zj = Weight{j} * A{j-1};
            Aj = sigmoid(Zj);
            A{j} = [ones(1, size(Aj, 2)); Aj];
        end
        Weight{nlayer+1} = reshape(Wlayed(W{nlayer+1, 2}(1):W{nlayer+1, 2}(2)), W{nlayer+1, 1});
        Znlayer1 = Weight{nlayer+1} * A{nlayer};
        A{nlayer+1} = sigmoid(Znlayer1);
        yforward = A{nlayer+1};
        
        % Computation of derivatives:
        error = cell(nlayer+1, 1);
        gradJ = cell(nlayer+1, 1);
        error{nlayer+1} = yforward - ytrain;
        Wreg = [zeros(size(Weight{nlayer+1}, 1), 1), Weight{nlayer+1}(:, 2:end)];
        gradJ{nlayer+1} = (1 / xm) * (error{nlayer+1} * A{nlayer}' + lambda * Wreg); 
        for j = fliplr(2:nlayer)
            errorj = (Weight{j+1}' * error{j+1}) .* A{j} .* (1 - A{j});
            error{j} = errorj(2:end, :);
            Wreg = [zeros(size(Weight{j}, 1), 1), Weight{j}(:, 2:end)];
            gradJ{j} = (1 / xm) * (error{j} * A{j-1}' + lambda * Wreg); 
        end
        error1 = (Weight{2}' * error{2}) .* A{1} .* (1 - A{1});
        error{1} = error1(2:end, :);
        Wreg = [zeros(size(Weight{1}, 1), 1), Weight{1}(:, 2:end)];
        gradJ{1} = (1 / xm) * (error{1} * Xtrain' + lambda * Wreg);
        
        % Total gradient computation:
        gradT = zeros(size(Wlayed));
        for j = 1:(nlayer+1)
            gradT(W{j, 2}(1):W{j, 2}(2)) = gradJ{j}(:);
        end
        
        % Update:
        Loss{1} = Loss{2};
        Loss{2} = gradT;
        
        % Sign of derivate multiplication:
        losstemp = Loss{1} .* Loss{2};
        
        % Delta & which has already reached maximum:
        Delta = mu_plus*Delta0.*(losstemp > 0) + mu_moins*Delta0.*(losstemp < 0) ...
                        + Delta0.*(losstemp == 0);
        Delta(Delta > deltamax) = deltamax;
        Delta0 = Delta;
        
        % Updating weights:
        DeltaW = -Delta0.*(Loss{2} > 0) + Delta0.*(Loss{2} < 0);
        Wlayed = Wlayed + DeltaW;
        
        % Info on the screen:
        if mod(i, epoch/10) == 0
            fprintf("Iteration %i out of %i. \n", [i, epoch]);
        end

    end
    
    % We return the weights:
    out = cell(nlayer+1, 1);
    for i = 1:(nlayer+1)
        out{i} = reshape(Wlayed(W{i, 2}(1):W{i, 2}(2)), W{i, 1});
    end

end



function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end