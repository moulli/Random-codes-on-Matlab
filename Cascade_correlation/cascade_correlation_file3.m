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

% Redefining data set
Xtrain = dataTempShuffle(1:trainSize, 1:(end-1))';
ytrain = dataTempShuffle(1:trainSize, end)';
Xtest = dataTempShuffle((trainSize+1):end, 1:(end-1))';
ytest = dataTempShuffle((trainSize+1):end, end)';



%% Testing neural net function

structure = [20 20];
lambda = 3;
epoch = 1000;
W = f_RpropNNTrain(Xtrain, ytrain, structure, lambda, epoch);

% Precision on training set:
yguess = f_RpropNNTest(W, Xtrain);
precision = sum(yguess == ytrain) / length(yguess);
fprintf("Precision on training set with neural net: %f \n", precision);

% Precision on test set:
yguess = f_RpropNNTest(W, Xtest);
precision = sum(yguess == ytest) / length(yguess);
fprintf("Precision on test set with neural net: %f \n \n", precision);



%% Testing cascade correlation algorithm

pool = 10;
time = 0;
sen = 10^(-6);
[Wnet, arch] = f_RpropCasCorTrain(Xtrain, ytrain, lambda, pool, time, sen);
[youttrain, prectrain] = f_RpropCasCorTest(Xtrain, ytrain, Wnet, arch);
[youttest, prectest] = f_RpropCasCorTest(Xtest, ytest, Wnet, arch);

fprintf("Precision on training set with CCNN: %f \n", prectrain);
fprintf("Precision on test set with CCNN: %f \n", prectest);



%% Testing cascade correlation with Polak-Ribière:

pool = 10;
time = 1;
sen = 10^(-1);
[Wnet, arch] = f_PolakCasCorTrain(Xtrain, ytrain, lambda, pool, time, sen);
[youttrain, prectrain] = f_PolakCasCorTest(Xtrain, ytrain, Wnet, arch);
[youttest, prectest] = f_PolakCasCorTest(Xtest, ytest, Wnet, arch);

fprintf("Precision on training set with CCNN: %f \n", prectrain);
fprintf("Precision on test set with CCNN: %f \n", prectest);



%% Functions

function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end