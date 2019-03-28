clear; close all; clc



%% Loading data

load monkeydata_training.mat

% Basic informations:
[ntest, nmov] = size(trial);
nneu = size(trial(1, 1).spikes, 1);



%% Definition of a training set

% Matrix with number of spikes per neuron, per recording:
Xspikes = zeros(ntest*nmov, nneu);
for j = 1:nmov
    for i = 1:ntest
        temp = trial(i, j).spikes(:, 1:300);
        nsp = sum(temp, 2);
        Xspikes(i + ntest*(j-1), :) = nsp';
    end
end

% Vector of movements associated to these spikes vectors:
yspikes = zeros(ntest*nmov, 1);
for i = 1:nmov
    temp = i * ones(ntest, 1);
    yspikes((ntest*(i-1)+1):(ntest*i)) = temp;
end

% Shuffling the data:
percentage = 0.80;
dataTemp = [Xspikes, yspikes];
shuffleVect = randperm(ntest*nmov);
dataTempShuffle = dataTemp(shuffleVect, :);
trainSize = percentage * ntest * nmov;
Xtrain = dataTempShuffle(1:trainSize, 1:(end-1))';
ytrain = dataTempShuffle(1:trainSize, end)';
Xtest = dataTempShuffle((trainSize+1):end, 1:(end-1))';
ytest = dataTempShuffle((trainSize+1):end, end)';


Xtrain = [ones(1, size(Xtrain, 2)); Xtrain];
Xtest = [ones(1, size(Xtest, 2)); Xtest];
ytrain_n = ytrain == [1:8]';
m = size(Xtrain, 2);


% Definition of initial weights:
W0 = 0.012 * (rand(nmov, nneu + 1)-0.5);

% Definition of parameters for resilient backprop:
mu_plus = 1.2;
mu_moins = 0.5;
deltamax = 50;
epoch = 40;
Delta0 = 0.01 * ones(size(W0));
Loss = {zeros(size(W0)), zeros(size(W0))};

% Optimization using resilient backpropagation:
for i = 1:epoch
    
    % Forward propagation:
    yforward = sigmoid(W0*Xtrain);
    % Computation of derivatives:
    error = yforward - ytrain_n;
    gradJ0 = (1 / m) * error * Xtrain';
    % Update:
    Loss{1} = Loss{2};
    Loss{2} = gradJ0;
    % Sign of derivate multiplication:
    losstemp = Loss{1} .* Loss{2};
    % Delta & which has already reached maximum:
    Delta = mu_plus*Delta0.*(losstemp > 0) + mu_moins*Delta0.*(losstemp < 0) ...
                    + Delta0.*(losstemp == 0);
    Delta(Delta > deltamax) = deltamax;
    Delta0 = Delta;
    % Updating weights:
    DeltaW0 = -Delta0.*(Loss{2} > 0) + Delta0.*(Loss{2} < 0);
    W0 = W0 + DeltaW0;
    
end
% Saving W0
W0save = W0;


% Testing precision on training set:
A = sigmoid(W0*Xtrain);
yguess = zeros(1, size(A, 2));
for i = 1:length(yguess)
    [mind, ind] = max(A(:, i));
    yguess(i) = ind;
end
precision = sum(yguess == ytrain) / length(yguess);
fprintf("Precision on training set: %f \n", precision);

% Testing precision on test set:
Atest = sigmoid(W0*Xtest);
yguess = zeros(1, size(Atest, 2));
for i = 1:length(yguess)
    [mind, ind] = max(Atest(:, i));
    yguess(i) = ind;
end
precision = sum(yguess == ytest) / length(yguess);
fprintf("Precision on test set: %f \n \n", precision);


neuron = 0.012 * (rand(1, nneu + 1)-0.5);
for i = 1:epoch
    
end






function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end