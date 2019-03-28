clear; close all; clc



%% Construction of the dataset:
x = 0:0.5:10;
y = 0:0.5:10;
trol = length(x);
[X, Y] = meshgrid(x, y);
t = exp(-X.^2 -(Y-5).^2 + X.*Y + X - 1);
t = int8(t > 300);
% figure
% subplot(1, 2, 1)
% surf(t)
Xdata = [X(:), Y(:)];
ydata = t(:);
% Check
X = reshape(Xdata(:, 1), trol, trol);
Y = reshape(Xdata(:, 2), trol, trol);
t = reshape(ydata, trol, trol);
% subplot(1, 2, 2)
% surf(t)

% NEW DATA:
trol = 4;
Xdata = [1, 1; 1, 2; 1, 3; 1, 4; 2, 1; 2, 2; 2, 3; 2, 4; 3, 1; 3, 2; 3, 3; 3, 4; 4, 1; 4, 2; 4, 3; 4, 4];
ydata = [1; 1; 1; 1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
t = reshape(ydata, trol, trol);


%% Normalizing data:
Xdata = (Xdata - mean(Xdata)) ./ std(Xdata);



%% Training & test set and constants:
prop = 0.8;
shuffle = randperm(size(Xdata, 1));
Xdata = Xdata(shuffle, :);
ydata = ydata(shuffle, :);
Xdata = Xdata';
ydata = ydata';
Xtrain = Xdata(:, 1:floor(prop*size(Xdata, 2)));
ytrain = ydata(1:floor(prop*size(Xdata, 2)));
Xtest = Xdata(:, floor(prop*size(Xdata, 2))+1:end);
ytest = ydata(floor(prop*size(Xdata, 2))+1:end);

n = size(Xdata, 1) - 1;
m = size(Xdata, 2);
class = length(unique(ydata));


%% Training with only the linear relationship:
% Definition of initial weights:
W0 = 0.012 * (rand(class, n + 1)-0.5);
yinit = sigmoid(W0*Xdata);



% Cascade correlation:
lambda = 0;
sen = 0.1;
pool = 10;

% W0 = f_resilientBackProp(Xtrain, ytrain, W0, lambda, sen);
% nbsib = 1;
% [theta, type] = f_maxCovariance(Xtrain, ytrain, W0, pool, sen, nbsib);

time = 3;
[Wnet, architecture] = f_RpropCasCorTrain(Xtrain, ytrain, lambda, pool, time, sen);
[yout, prec] = f_RpropCasCorTest(Xdata, ydata, Wnet, architecture);
prec

% Classic resilient backprop:
epoch = 10;
Wrb = f_RpropNNTrain(Xtrain, ytrain, [10], lambda, epoch);
yrb = f_RpropNNTest(Wrb, Xdata);
prec_rb = sum(yrb == ydata) / length(ydata);
fprintf("Precision with regular NN and resilient backpropagation: %f \n", prec_rb);





% % Testing precision on training set:
% A = sigmoid(W0*Xtrain);
% [ma, yguess] = max(A);
% yguess = yguess-1;
% [yguess', ytrain'];
% precision = sum(yguess == ytrain) / length(yguess);
% fprintf("Precision on training set: %f \n", precision);
% 
% % Testing precision on test set:
% Atest = sigmoid(W0*Xtest);
% [ma, yguess] = max(Atest);
% yguess = yguess-1;
% precision = sum(yguess == ytest) / length(yguess);
% fprintf("Precision on test set: %f \n \n", precision);

% Plotting predictions:
Adatainit = yinit;
[mainit, ydatainit] = max(Adatainit);
ydatainit = ydatainit-1;
DATAinit = [Xdata', ydatainit'];
Xtempinit(shuffle, :) = DATAinit;
Xtempinit = reshape(Xtempinit(:, 3), trol, trol);

% Adata = sigmoid(W0*Xdata);
% [ma, ydata] = max(Adata);
% ydata = ydata-1;
% DATA = [Xdata', ydata'];
% Xtemp(shuffle, :) = DATA;
% Xtemp = reshape(Xtemp(:, 3), trol, trol);

DATAout = [Xdata', yout'];
Xtempout(shuffle, :) = DATAout;
Xtempout = reshape(Xtempout(:, 3), trol, trol);

DATArb = [Xdata', yrb'];
Xrb(shuffle, :) = DATArb;
Xrb = reshape(Xrb(:, 3), trol, trol);


figure
subplot(1, 4, 1)
surf(t)
title("Actual data")
subplot(1, 4, 2)
surf(Xtempinit)
title("Guess before training")
% subplot(1, 4, 3)
% surf(Xtemp)
subplot(1, 4, 3)
surf(Xtempout)
title("Guess after training")
subplot(1, 4, 4)
surf(Xrb)
title("Guess with resilient backprop")


% i = 0;
% tic
% while toc < 1
%     i = i + 1;
% end
% toc, i



        










function g = sigmoid(z)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

    g = 1.0 ./ (1.0 + exp(-z));
end