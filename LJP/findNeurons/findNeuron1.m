clear; close all; clc
addpath('C:/Users/Hippolyte Moulle/Desktop/Amusements_Matlab/Heavy_files')


load('mstd.mat')

%% Plotting matrix:
couche = 5;
figure
image(mstd(:, :, couche), 'CDataMapping', 'scaled')
colorbar
axis equal

%% No noise:
figure
image(mstd(:, :, couche) .* (mstd(:, :, couche) > 15), 'CDataMapping', 'scaled')
colorbar
axis equal

%% Isolating neurons:
yneu = [566, 583, 572, 612, 848, 788, 388, 782];
xneu = [180, 415, 436, 423, 517, 474, 299, 352];
sxneu = [5, 6, 5, 7, 7, 6, 6, 6];
syneu = [5, 6, 7, 8, 8, 8, 9, 5];
% figure
% for i = 1:8
%     subplot(4, 2, i)
%     image(mstd((xneu(i)-ceil(sxneu(i)/2)):(xneu(i)+ceil(sxneu(i)/2)), (yneu(i)-ceil(syneu(i)/2)):(yneu(i)+ceil(syneu(i)/2)), couche), 'CDataMapping', 'scaled')
%     colorbar
%     axis equal
% end

%% Creation of fake neuron:
meaneu = [0, 0];
covneu = [7, -3; -3, 7];
maxneu = 120 * sqrt((2*pi)^2*det(covneu));
coord = -3:3;
[X, Y] = meshgrid(coord, coord);
neureq = maxneu * mvnpdf([X(:), Y(:)], meaneu, covneu);
% figure
% image(reshape(neureq, size(X)), 'CDataMapping', 'scaled')
% colorbar

%% Creating whole plan of fake neurons:
meanx = 4.5; stdx = 1;
meany = 4.5; stdy = 1;
covcor = 0; stdcov = 0.5;
meanoise = 5; stdnoise = 5;
planpar = [250, 100, 50, 15;
           50, 30, 20, 4;
           0, 0.1, 1.5, 2.5];
splan = size(planpar, 2);
numneur = 10000; numstd = 1000;
numfin = round(numneur + numstd*randn);
coordfin = zeros(numfin, 2);
plan = zeros([size(mstd(:, :, 1), 1), size(mstd(:, :, 1), 2), splan]);

[Xfield, Yfield] = meshgrid(-6:6, -6:6);
Xfield = Xfield(:);
Yfield = Yfield(:);
meaneu = [0, 0];
for i = 1:numfin
    planrand = splan*rand;
    planchoice = find(planrand <= [planpar(3, :), splan]);
    planchoice = planchoice(1) - 1;
    finding = -1;
    numround = 0;
    while finding == -1 && numround < 100
        randx = round(rand * (size(mstd, 1)-15))+7;
        randy = round(rand * (size(mstd, 2)-15))+7;
        anafield = plan((randx-6):(randx+6), (randy-6):(randy+6), planchoice);
        if nnz(anafield) == 0
            break
        end
        numround = numround + 1;
    end   
    if numround == 100
        continue
    end
    numround = 0;
    matcorval = covcor + stdcov*randn;
    matrixstd = [meanx + stdx*randn, matcorval;
                 matcorval, meany + stdy*randn];
    if det(matrixstd) < 0
        continue
    else
        coordfin(i, 1) = randx;
        coordfin(i, 2) = randy;
    end
    numval = (planpar(1, planchoice) + planpar(2, planchoice)*randn) * sqrt((2*pi)^2*det(matrixstd));
    valexp = numval * mvnpdf([Xfield, Yfield], meaneu, matrixstd);
    plan((randx-6):(randx+6), (randy-6):(randy+6), planchoice) = reshape(valexp, [13, 13]);
end

plantemp = zeros(numel(mstd(:, :, 1)), splan+1);
noise = meanoise + stdnoise*randn(numel(mstd(:, :, 1)), 1);
plantemp(:, end) = noise;
for i = 1:splan
    planttemp = plan(:, :, i);
    plantemp(:, i) = planttemp(:);
end
planfinal = reshape(max(plantemp, [], 2), size(mstd(:, :, 1)));
coordfin = coordfin(all(coordfin, 2), :);
figure
image(planfinal, 'CDataMapping', 'scaled')
colorbar
axis equal






