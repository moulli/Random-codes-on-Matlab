clear; close all; clc




% Data implementation
RGB255 = imread('champ.jpg');
figure
subplot(1, 2, 1)
image(RGB255)
title("Initial image", "Interpreter", "latex")
RGB = rgb2ntsc(RGB255);
% Weird way of rearranging the colors, but usable matrix
[h, l, c] = size(RGB);
Red = RGB(:, :, 1);
Green = RGB(:, :, 2);
Blue = RGB(:, :, 3);
data0 = [Red(:), Green(:), Blue(:)];
t = size(data0, 1);

% Centroids initialization
clus = 6;
centroid = zeros(clus, c);
for i = 1:clus
    n = floor(rand(1)*t);
    centroid(i, :) = data0(n, :);
end
% Assignment of the points to the nearest centroid:
for i = 1:t
    D = sqrt(sum((data0(i, :) - centroid).^2, 2));
    [M, index] = min(D);
    Ind(i) = index;
end

% K-means algorithm
num_iter = 15;
for i = 1:num_iter
    for j = 1:clus
        top = (Ind == j);
        ntop = sum(top);
        if ntop == 0
            centroid(j, :) = [10000, 10000];
        else
            centroid(j, :) = mean(data0(top, :));
        end
    end
    for j = 1:t
        D = sqrt(sum((data0(j, :) - centroid).^2, 2));
        [M, index] = min(D);
        Ind(j) = index;
    end
    fprintf("Iteration %d \n", i);
end

data1 = zeros(size(data0));
for i = 1:clus
    data1(Ind == i, :) = (centroid(i, :) .* ones(size(data1(Ind == i, :))));
end
Red1 = reshape(data1(:, 1), h, l);
Green1 = reshape(data1(:, 2), h, l);
Blue1 = reshape(data1(:, 3), h, l);
RGB1 = zeros(h, l, c);
    RGB1(:, :, 1) = Red1;
    RGB1(:, :, 2) = Green1;
    RGB1(:, :, 3) = Blue1;

% Displaying the new image
RGB1255 = ntsc2rgb(RGB1);
subplot(1, 2, 2)
image(RGB1255)
title("Image simplified, using centroids", "Interpreter", "latex")

