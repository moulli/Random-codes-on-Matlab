clear; close all; clc



x = -1.5:0.01:1.5;
y = -0.5:0.01:1.5;
[X, Y] = meshgrid(x, y);
r = (X - 1).^2 + 10 * (X.^2 - Y).^2;

K = 8;
figure
subplot(2, K, 1)
contour(x, y, r, [0.01 0.1 0.4 1:20]);
grid on
hold on

% Control of convergence:
power = 5;

% First shot with the points:
p = 25;
xrand = rand(2*p, 1)*(x(end) - x(1)) + x(1);
yrand = rand(2*p, 1)*(y(end) - y(1)) + y(1);
scatter(xrand, yrand, "+k")
rrand = (xrand - 1).^2 + 10 * (xrand.^2 - yrand).^2;
D = 1 ./ rrand.^power;
% % 1st way, completely fucked up:
% Dp = round(D ./ min(D));
% xtemp = zeros(sum(Dp), 2);
% for i = 1:length(Dp)
%     xtemp(1+sum(Dp(1:i-1)):sum(Dp(1:i)), :) = ones(Dp(i), 2) .* [xrand(i), yrand(i)];
% end
% mu = mean(xtemp)
% sigma = cov(xtemp)
% 2nd way, more mathematical:
Dm = D ./ sum(D);
Xrand = [xrand, yrand];
mu = sum(Xrand .* Dm);
sigma = (Xrand - mu)' * diag(Dm) * (Xrand - mu);


% Following shots with the points:
for i = 2:K
    Prand1 = mvnrnd(mu, sigma, p);
    xrand1 = Prand1(:, 1);
    yrand1 = Prand1(:, 2);
    xrand2 = rand(p, 1)*(x(end) - x(1)) + x(1);
    yrand2 = rand(p, 1)*(y(end) - y(1)) + y(1);
    xrand = [xrand1; xrand2];
    yrand = [yrand1; yrand2];
    subplot(2, K, i)
    contour(x, y, r, [0.01 0.1 0.4 1:20]);
    grid on
    hold on
    scatter(xrand1, yrand1, "+r")
    hold on
    scatter(xrand2, yrand2, "+k")
    rrand = (xrand - 1).^2 + 10 * (xrand.^2 - yrand).^2;
    D = 1 ./ rrand.^power;
    Dm = D ./ sum(D);
    Xrand = [xrand, yrand];
    mu = sum(Xrand .* Dm);
    sigma = (Xrand - mu)' * diag(Dm) * (Xrand - mu);
end

% Continuing:
for i = 1:K-1
    Prand1 = mvnrnd(mu, sigma, p);
    xrand1 = Prand1(:, 1);
    yrand1 = Prand1(:, 2);
    xrand2 = rand(p, 1)*(x(end) - x(1)) + x(1);
    yrand2 = rand(p, 1)*(y(end) - y(1)) + y(1);
    xrand = [xrand1; xrand2];
    yrand = [yrand1; yrand2];
    subplot(2, K, K+i)
    contour(x, y, r, [0.01 0.1 0.4 1:20]);
    grid on
    hold on
    scatter(xrand1, yrand1, "+r")
    hold on
    scatter(xrand2, yrand2, "+k")
    rrand = (xrand - 1).^2 + 10 * (xrand.^2 - yrand).^2;
    D = 1 ./ rrand.^power;
    Dm = D ./ sum(D);
    Xrand = [xrand, yrand];
    mu = sum(Xrand .* Dm);
    sigma = (Xrand - mu)' * diag(Dm) * (Xrand - mu);
end

% Last point (we don't want to have random points bias the final
% coordinates):
Prand1 = mvnrnd(mu, sigma, 2*p);
xrand = Prand1(:, 1);
yrand = Prand1(:, 2);
subplot(2, K, 2*K)
contour(x, y, r, [0.01 0.1 0.4 1:20]);
grid on
hold on
scatter(xrand, yrand, "+r")
rrand = (xrand - 1).^2 + 10 * (xrand.^2 - yrand).^2;
D = 1 ./ rrand.^power;
Dm = D ./ sum(D);
Xrand = [xrand, yrand];
mu = sum(Xrand .* Dm);
sigma = (Xrand - mu)' * diag(Dm) * (Xrand - mu);


fprintf("Minimum of function is evaluated at coordinates x: %f & y: %f \n", ...
            [mu(1), mu(2)]);



