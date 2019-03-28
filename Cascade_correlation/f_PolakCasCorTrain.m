function [Wnet, architecture] = f_PolakCasCorTrain(Xtrain, ytrain, lambda, pool, time, sen)
% Train neural network based on cascade correlation algorithm.
% Xtrain is the data set, ytrain labels.
% lambda is the regularization parameter.
% pool is the number of siblings&descendents we want to try at each iteration.
% time is the time we want the algorithm to train.
% sen is the sensibility at which we want to skip to next task:
% (eg add a neuron or stop training weights).
% W is a cell is which the weights are located.
% architecture is a vector that indicates to which layer each neuron belongs.
% Hippolyte MOULLE


    %% First indication for time:
    tic


    %%  Initialization:
    Xtrain = [ones(1, size(Xtrain, 2)); Xtrain];
    % Reconstruction of a usable label vector:
    if size(ytrain, 1) > 1
        error("Output vector must be provided as a row vector")
    end
    class = unique(ytrain);
    yreach = (ytrain == (class(1):class(end))');
    
    
    %% Initial perceptron:
    n = size(Xtrain, 1);
    label = length(class(1):class(end));
    winit = 0.012 * (rand(label, n)-0.5);
    W = f_PolakResilientBasic(Xtrain, ytrain, winit, lambda, sen);
    % Initial network weights and architecture:
    Wnet = {W};
    architecture = [0];
    
    
    %% Adding first neuron:
    % Finding best neuron to add to the network (0 siblings already):
    [theta, type] = f_PolakMaxCovariance(Xtrain, ytrain, W, pool, sen, 0);
    % Freezing its weights and its location in the net architecture:
    Wnet{end+1} = theta;
    architecture(end+1) = 1; % it has to be a descendent
    % Building a new exploitable training set:
    Xtrain = [Xtrain; sigmoid(theta * Xtrain)];
    % Modifying weight matrix to include new column:
    W = [W, 0.012 * (rand(label, 1)-0.5)];
%     W = 0.012 * (rand(size(W, 1), size(W, 2) + 1)-0.5);
    % Training the newly added weights:
    W = f_PolakResilientBasic(Xtrain, ytrain, W, lambda, sen);
    
    
    %% Adding neurons one after the other and training the weights:
    tic % start time
    nbsib = 1; % number of siblings already in the layer (1 by default)
    while toc < time
        % Finding best neuron to add to the network:
        [theta, type] = f_PolakMaxCovariance(Xtrain, ytrain, W, pool, sen, nbsib);
        % Freezing its weights and its location in the net architecture:
        Wnet{end+1} = theta;
        architecture(end+1) = architecture(end) + type;
        % Modifyin number of siblings:
        if type == 0 
            nbsib = nbsib + 1;
        elseif type == 1
            nbsib = 1;
        end
        % Building a new exploitable training set:
        Xtrain = [Xtrain; sigmoid(theta * Xtrain(1:end-nbsib+1, :))];
        % Modifying weight matrix to inclue new row:
%         W = [W, 0.012 * (rand(label, 1)-0.5)];
        W = 0.012 * (rand(size(W, 1), size(W, 2) + 1)-0.5);
        % Training the newly added weights:
        W = f_PolakResilientBasic(Xtrain, ytrain, W, lambda, sen);
    end
    
    
    %% Updating first matrix in Wnet:
    Wnet{1} = W;
    % ATTENTION: all the unfreezed weights are stocked in Wnet{1}!


end



function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end