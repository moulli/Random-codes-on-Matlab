clear; close all; clc

%% Cell life algorithm:
% Based on an initial environment, with nline lines and ncol columns.
% Cells are randomly initialized.
% If a cell is surrounded by less than mincell, it dies.
% If a cell is surrounded by more than maxcell, it dies.
% Otherwise, a new cell is created, randomly, near the old cell.


%% Parameters:
nline = 20;
ncol = 20;
epoch = 2;
timer = 0.000001;
mincell = 2;
maxcell = 3;


%% Initializing environment:
envi = round(rand(nline, ncol)); % random initialization
envi = zeros(nline, ncol);
envi(floor(nline/3)-floor(nline/10):floor(nline/3)+floor(nline/9), ...
     floor(ncol/3)-floor(ncol/5):floor(ncol/3)+floor(ncol/8)) = 1; %location
envi(floor(nline/2)+floor(nline/4), floor(ncol/2)-1:floor(ncol/2)+1) = 1;
envi(floor(nline/2)-1+floor(nline/4), floor(ncol/2)+1) = 1;
envi(floor(nline/2)-2+floor(nline/4), floor(ncol/2)) = 1;
enviI = envi;


%% Algorithm:
figure
image(envi, "CDataMapping", "scaled")
title('Cell population, initial')
blackwhite = [(0.8:-0.001:0.2)', (0.8:-0.001:0.2)', (0.8:-0.001:0.2)']; 
colormap(gca, blackwhite);
c = colormap;

for k = 1:epoch
    % Random update:
    R = randperm(nline*ncol);
    for d = 1:nline*ncol
        i = mod(R(d), nline);
        if i == 0
            i = nline;
        end
        j = ceil(R(d) / nline);
        if sum(envi(:)) == 0
            break
        end
        [indice, value] = fNeighbours(envi, i, j);
        if sum(value) < mincell || sum(value) > maxcell
            envi(i, j) = 0;
        elseif sum(value) == 3
            envi(i, j) = 1;
        end
    end
    pause(timer);
    image(envi, "CDataMapping", "scaled")
    title(['Cell population, epoch: ', num2str(k)])
    blackwhite = [(0.8:-0.001:0.2)', (0.8:-0.001:0.2)', (0.8:-0.001:0.2)']; 
    colormap(gca, blackwhite);
    c = colormap;
end




% for k = 1:epoch
%     % Random update:
%     R = randperm(nline*ncol);
%     for d = 1:nline*ncol
%         i = mod(R(d), nline);
%         if i == 0
%             i = nline;
%         end
%         j = ceil(R(d) / nline);
%         if sum(envi(:)) == 0
%             break
%         end
%         if envi(i, j) == 1
%             [indice, value] = fNeighbours(envi, i, j);
%             if sum(value) >= mincell && sum(value) <= maxcell
%                 nullval = (value == 0) .* 1;
%                 potential = indice .* nullval;
%                 onepot = find(potential == 1);
%                 if isempty(onepot) == 1
%                     break
%                 end
%                 select = randperm(length(onepot));
%                 keep = onepot(select(1));
%                 [it, jt] = fChange(keep);
%                 envi(i+it, j+jt) = 1;
%             else
%                 envi(i, j) = 0;
%             end
%         end
%     end
%     pause(timer);
%     image(envi, "CDataMapping", "scaled")
%     title(['Cell population, epoch: ', num2str(k)])
%     blackwhite = [(0.8:-0.001:0.2)', (0.8:-0.001:0.2)', (0.8:-0.001:0.2)']; 
%     colormap(gca, blackwhite);
%     c = colormap;
% end



% for k = 1:epoch
%     % Random update:
%     iR = randperm(nline);
%     jR = randperm(ncol);
%     % Death of cells:
%     enviD = zeros(nline, ncol);
%     for i = 1:nline
%         for j = 1:ncol
%             if envi(iR(i), jR(j)) == 1
%                 [~, value] = fNeighbours(envi, iR(i), jR(j));
%                 if sum(value) >= mincell && sum(value) <= maxcell
%                     enviD(iR(i), jR(j)) = 1;
%                 end
%             end
%         end
%     end
%     % Creation of new cells:
%     enviN = enviD;
%     for i = 1:nline
%         for j = 1:ncol
%             if enviD(iR(i), jR(j)) == 1
%                 [indice, value] = fNeighbours(enviN, iR(i), jR(j));
%                 nullval = (value == 0) .* 1;
%                 potential = indice .* nullval;
%                 onepot = find(potential == 1);
%                 if isempty(onepot) == 1
%                     break
%                 end
%                 select = randperm(length(onepot));
%                 keep = onepot(select(1));
%                 [it, jt] = fChange(keep);
%                 enviN(iR(i)+it, jR(j)+jt) = 1;
%             end
%         end
%     end
%     envi = enviN;
%     pause(timer);
%     image(envi, "CDataMapping", "scaled")
%     title(['Cell population, epoch: ', num2str(k)])
%     blackwhite = [(0.8:-0.001:0.2)', (0.8:-0.001:0.2)', (0.8:-0.001:0.2)']; 
%     colormap(gca, blackwhite);
%     c = colormap;
% end
                    




% for k = 1:epoch
%     enviN = zeros(nline+2, ncol+2);
%     enviT = zeros(nline+2, ncol+2);
%     enviT(2:end-1, 2:end-1) = envi;
%     enviT(2:end-1, 1) = envi(:, end);
%     enviT(2:end-1, end) = envi(:, 1);
%     enviT(1, :) = enviT(end-1, :);
%     enviT(end, :) = enviT(2, :);
%     for i = 2:nline+1
%         for j = 2:ncol+1
%             entou = enviT(i-1:i+1, j-1:j+1);
%             if sum(entou(:)) > 1 && sum(entou(:)) < 6
%                 % Keeping the cell:
%                 enviN(i, j) = 1;
%                 % Creating new random cell:
%                 newcell = 0;
%                 count = 0;
%                 while newcell == 0 && count < 100
%                     inew = randperm(3)-2;
%                     it = inew(1);
%                     jnew = randperm(3)-2;
%                     jt = jnew(1);
%                     if enviT(i+it, j+jt) == 0 && enviN(i+it, j+jt) == 0
%                         newcell = 1;
%                     end
%                     count = count + 1;
%                 end
%                 enviN(i+it, j+jt) = 1;
%                 enviN(2, :) = enviN(end, :);
%                 enviN(end-1, :) = enviN(1, :);
%                 enviN(:, 2) = enviN(:, end);
%                 enviN(:, end-1) = enviN(:, 1);
%             end
%         end
%     end
%     envi = enviN(2:end-1, 2:end-1);
%     pause(timer);
%     image(envi, "CDataMapping", "scaled")
%     title(['Cell population, epoch: ', num2str(k)])
%     blackwhite = [(0.8:-0.001:0.2)', (0.8:-0.001:0.2)', (0.8:-0.001:0.2)']; 
%     colormap(gca, blackwhite);
%     c = colormap;
% end

figure
subplot(1, 2, 1)
image(enviI, "CDataMapping", "scaled")
title('Cell population, epoch: 0')
blackwhite = [(0.8:-0.001:0.2)', (0.8:-0.001:0.2)', (0.8:-0.001:0.2)']; 
colormap(gca, blackwhite);
c = colormap;
subplot(1, 2, 2)
image(envi, "CDataMapping", "scaled")
title('Cell population, epoch: end')
blackwhite = [(0.8:-0.001:0.2)', (0.8:-0.001:0.2)', (0.8:-0.001:0.2)']; 
colormap(gca, blackwhite);
c = colormap;
    







function [indice, value] = fNeighbours(A, i, j)
% Function that returns the neighbours of a value in a matrix.
% indice is under the form:
% 1 | 2 | 3
% 4 | * | 5
% 6 | 7 | 8
% or:
% * | 5
% 7 | 8
% or:
% 1 | 2 
% 4 | *
% And value is the values associated, for instance [1, 1, 0, 1, 0, 0, 0, 0].
    [n, m] = size(A);
    if i == 1 && j == 1
        indice = [0, 0 ,0, 0, 1, 0, 1, 1];
        value = [0, 0, 0, 0, A(i, j+1),  0, A(i+1, j), A(i+1, j+1)];
    elseif i == 1 && j == m
        indice = [0, 0 ,0, 1, 0, 1, 1, 0];
        value = [0, 0, 0, A(i, j-1), 0,  A(i+1, j-1), A(i+1, j), 0];
    elseif i == n && j == 1
        indice = [0, 1 ,1, 0, 1, 0, 0, 0];
        value = [0, A(i-1, j), A(i-1, j+1), 0, A(i, j+1),  0, 0, 0];
    elseif i == n && j == m
        indice = [1, 1 ,0, 1, 0, 0, 0, 0];
        value = [A(i-1, j-1), A(i-1, j), 0, A(i, j-1), 0,  0, 0, 0];
    elseif i == 1 && j ~= 1 && j ~= m
        indice = [0, 0 ,0, 1, 1, 1, 1, 1];
        value = [0, 0, 0, A(i, j-1), A(i, j+1),  ...
                 A(i+1, j-1), A(i+1, j), A(i+1, j+1)];
    elseif i == n && j ~= 1 && j ~= m
        indice = [1, 1 ,1, 1, 1, 0, 0, 0];
        value = [A(i-1, j-1), A(i-1, j), A(i-1, j+1), A(i, j-1), A(i, j+1), ...
                 0, 0, 0];
    elseif i ~= 1 && i ~= n && j == 1
        indice = [0, 1 ,1, 0, 1, 0, 1, 1];
        value = [0, A(i-1, j), A(i-1, j+1), 0, A(i, j+1),  ...
                 0, A(i+1, j), A(i+1, j+1)];
    elseif i ~= 1 && i ~= n && j == m
        indice = [1, 1 ,0, 1, 0, 1, 1, 0];
        value = [A(i-1, j-1), A(i-1, j), 0, A(i, j-1), 0,  ...
                 A(i+1, j-1), A(i+1, j), 0];
    else 
        indice = [1, 1, 1, 1, 1, 1, 1, 1];
        value = [A(i-1, j-1), A(i-1, j), A(i-1, j+1), ...
                 A(i, j-1), A(i, j+1), ...
                 A(i+1, j-1), A(i+1, j), A(i+1, j+1)];
    end
end

function [it, jt] = fChange(n)
% Function that returns the position of the new value based on n.
    if n == 1
        it = -1; jt = -1;
    elseif n == 2
        it = -1; jt = 0;
    elseif n == 3
        it = -1; jt = +1;
    elseif n == 4
        it = 0; jt = -1;
    elseif n == 5
        it = 0; jt = +1;
    elseif n == 6
        it = +1; jt = -1;
    elseif n == 7
        it = +1; jt = 0;
    elseif n == 8
        it = +1; jt = +1;
    end
end

