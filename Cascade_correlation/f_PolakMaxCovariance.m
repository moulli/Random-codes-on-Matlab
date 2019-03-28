function [theta, type] = f_PolakMaxCovariance(Xnew, y, Wnew, pool, sen, nbsib)
%% Function that performs maximum of covariance on a pool of neurons.
% Xnew is the set (using neurons already produced), y the labels.
% Wnew are the weights for perceptron already computed, and each added
% neuron's weight, so that it fits the Xnew set
% pool: number of neurons in each pool that will be trained independantly.
% sen is the sensibility of the error that will stop the algorithm.
% nbsib is tricky: it represents the number of siblings already in the layer.
% It will allow the program to not take into account all these neurons.
% By default, the user enters 1, because a layer was created just before.
% If 0 is entered, then it is the 1st neuron of the net, therefore this
% neuron is a descendent, and type must be 1.
% theta: weights; type: if neuron is sibling (0) or descendent (1).
% Hippolyte MOULLE 


    %%  Initialization:
    % Reconstruction of a usable label vector:
    if size(y, 1) > 1
        error("Output vector must be provided as a row vector")
    end
    class = unique(y);
    yreach = (y == (class(1):class(end))');
    % Pooling for sibling and descendent neurons:
    features = size(Xnew, 1);
    sibl = 0.012 * (rand(pool, features-nbsib)-0.5);
    desc = 0.012 * (rand(pool, features)-0.5);
    % Stocking the covariances:
    sibl_cov = zeros(pool, 1);
    desc_cov = zeros(pool, 1);


    %% Raising errors if dimension or sign problem:
    if size(Xnew, 2) ~= size(yreach, 2)
        error("Not same number of patterns in input and output")
    elseif size(Wnew, 1) ~= size(yreach, 1)
        error("Not same number of weights and classes")
    elseif size(Xnew, 1) ~= size(Wnew, 2)
        error("Different number of features for set and weights")
    elseif sen < 0
        error("Sensibility parameter must be positive")
    elseif pool < 0 || floor(pool) ~= pool
        error("Regularizaton parameter must be integer")
    end
    
    
    %% Error to be used to compute maximum of covariance:
    output = sigmoid(Wnew * Xnew);
    errare = output - yreach;
    errare = errare - mean(errare, 2);

    
    %% Parameters for while loop:
    % Definition of parameters for resilient backprop:
    mu_plus = 1.2;
    mu_moins = 0.5;
    deltamax = 50;
    % Exit function parameters:
    cov1 = 0;
    cov2 = 0;
    startCount = 0;
    
    
    %% Main algorithm:
    % We are going to work on each neuron in the pool one after the other.
    
    %% Sibling:
    for i = 1:pool
        % Updating resilient backprop parameters:
        Delta0 = 0.01 * ones(1, features-nbsib);
        Loss{1} = zeros(1, features-nbsib);
        Loss{2} = zeros(1, features-nbsib);
        while cov2 - cov1 > sen || startCount < 2
            % Computing neuron output:
            sibli = sibl(i, :);
            Vp = sigmoid(sibli * Xnew(1:end-nbsib, :));
            V = Vp - mean(Vp);
            % Updating covariance:
            cov1 = cov2;
            cov2 = sum(abs(sum(V .* errare, 2)));
            % 3d matrix to get rid of for loop:
            E3d = ones(size(errare, 1), size(errare, 2), features-nbsib) .* errare;
            % Computing gradient of covariance:
            dStemp = Xnew(1:end-nbsib, :) .* sigmoid(sibli * Xnew(1:end-nbsib, :)) ...
                        .* (1 - sigmoid(sibli * Xnew(1:end-nbsib, :)));
            dStemp = permute(dStemp, [3, 2, 1]); % to use 3d matrix
            dStemp = sum(abs(sum(dStemp .* E3d, 2)), 1);
            dSdw = squeeze(dStemp)'; % rearrange in a simple vector
            % Resilient backpropagation (MAXIMIZATION):
            % Update:
            Loss{1} = Loss{2};
            Loss{2} = dSdw;
            % Sign of derivate multiplication:
            losstemp = Loss{1} .* Loss{2};
            % Delta & which has already reached maximum:
            Delta = mu_plus*Delta0.*(losstemp > 0) + mu_moins*Delta0.*(losstemp < 0) ...
                            + Delta0.*(losstemp == 0);
            Delta(Delta > deltamax) = deltamax;
            Delta0 = Delta;
            % Updating weights (with a minus to maximize):
            DeltaW = -Delta0.*(Loss{2} > 0) + Delta0.*(Loss{2} < 0);
            sibl(i, :) = sibli - DeltaW;
            % Loop exit scheme:
            startCount = startCount + 1;
        end
        % Stocking retained covariance:
        Vp = sigmoid(sibl(i, :) * Xnew(1:end-nbsib, :));
        V = Vp - mean(Vp);
        sibl_cov(i) = sum(abs(sum(V .* errare, 2)));
        % Setting initial values again:
        cov1 = 0;
        cov2 = 0;
        startCount = 0;
    end
    
    %% Descendent:
    for i = 1:pool
        % Updating resilient backprop parameters:
        Delta0 = 0.01 * ones(1, features);
        Loss{1} = zeros(1, features);
        Loss{2} = zeros(1, features);
        while cov2 - cov1 > sen || startCount < 2
            % Computing neuron output:
            desci = desc(i, :);
            Vp = sigmoid(desci * Xnew);
            V = Vp - mean(Vp);
            % Updating covariance:
            cov1 = cov2;
            cov2 = sum(abs(sum(V .* errare, 2)));
            % 3d matrix to get rid of for loop:
            E3d = ones(size(errare, 1), size(errare, 2), features) .* errare;
            % Computing gradient of covariance:
            dStemp = Xnew .* sigmoid(desci * Xnew) .* (1 - sigmoid(desci * Xnew));
            dStemp = permute(dStemp, [3, 2, 1]); % to use 3d matrix
            dStemp = sum(abs(sum(dStemp .* E3d, 2)), 1);
            dSdw = squeeze(dStemp)'; % rearrange in a simple vector
            % Resilient backpropagation (MAXIMIZATION):
            % Update:
            Loss{1} = Loss{2};
            Loss{2} = dSdw;
            % Sign of derivate multiplication:
            losstemp = Loss{1} .* Loss{2};
            % Delta & which has already reached maximum:
            Delta = mu_plus*Delta0.*(losstemp > 0) + mu_moins*Delta0.*(losstemp < 0) ...
                            + Delta0.*(losstemp == 0);
            Delta(Delta > deltamax) = deltamax;
            Delta0 = Delta;
            % Updating weights (with a minus to maximize):
            DeltaW = -Delta0.*(Loss{2} > 0) + Delta0.*(Loss{2} < 0);
            desc(i, :) = desci - DeltaW;
            % Loop exit scheme:
            startCount = startCount + 1;
        end
        % Stocking retained covariance:
        Vp = sigmoid(desc(i, :) * Xnew);
        V = Vp - mean(Vp);
        desc_cov(i) = sum(abs(sum(V .* errare, 2)));
        % Setting initial values again:
        cov1 = 0;
        cov2 = 0;
        startCount = 0;
    end
    
    
    %% We are going to select the best neuron for the network:
    [ma, ind] = max([sibl_cov; desc_cov]);
    if ind <= pool
        theta = sibl(ind, :);
        type = 0;
    else
        theta = desc(ind - pool, :);
        type = 1;
    end


end



%% Annex sigmoid function:
function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end