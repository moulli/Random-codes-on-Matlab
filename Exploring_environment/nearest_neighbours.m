clear; close all; clc







%% Parameters

w = 10;
h = 10;
ax = [0, w, 0, h];
n = 200;
mnn = 30;
po = [1, 9];
pf = [9, 9];



%% Definition of the obstacles

lines = {};
polygons = {[3, 3; 3, 10;
             3, 10; 4, 10;
             4, 10; 4, 2;
             4, 2; 2, 2;
             2, 2; 2, 3;
             2, 3; 3, 3];
            [5.5, 0; 5.5, 8;
             5.5, 8; 6, 8;
             6, 8; 6, 0;
             6, 0; 5.5, 0];
            [7.5, 3; 7.5, 10;
             7.5, 10; 8.5, 10;
             8.5, 10; 8.5, 3;
             8.5, 3; 7.5, 3]};

         
         
%% Plot 

figure
for i = 1:size(lines, 1)
    p1 = lines{i, 1};
    p2 = lines{i, 2};
    line([p1(1), p2(1)], [p1(2), p2(2)], "Color", "black", "Linewidth", 2)
    hold on
end
for i = 1:size(polygons, 1)
    p1 = polygons{i};
    taille = size(p1, 1);
    p = [p1; p1(1, :)];
    for j = 1:taille
        line([p(j, 1), p(j+1, 1)], [p(j, 2), p(j+1, 2)], "Color", "black", "Linewidth", 2)
        hold on
    end
end
pt0 = rand(n-2, 2) .* [w, h];
pt0 = [po; pt0; pf];
pts = [];
for i = 1:n
    obs = 0;
    for j = 1:size(polygons, 1)
        p = polygons{j};
        if inpolygon(pt0(i, 1), pt0(i, 2), p(:, 1), p(:, 2)) == 1
            obs = 1;
        end
    end
    if obs == 0
        pts = [pts; pt0(i, :)];
    end
end
n = size(pts, 1);
scatter([po(1), pf(1)], [po(2), pf(2)], 10, [1, 0, 0])
hold on
scatter(pts(:, 1), pts(:, 2), 3)
title("Probabilistic Road Map, max. neighbours: 30", "Interpreter", "latex")
axis(ax)
hold on



%% Distance matrix

D = zeros(n);
for i = 1:n
    dist = sqrt(sum((pts - pts(i, :)) .^2, 2));
    D(:, i) = dist;
end
% Changing 0 to high value so that it does not influence distances
S = sum(D(:));
for i = 1:n
    D(i, i) = S;
end



%% Implementation of the PRM

% We have changed n, so we need to adapt mnn:
if mnn > n-1
    mnn = n-1;
end

% Computing the lines, avoiding obstacles
[mi, in] = sort(D);
in = in(1:(end-1), :); %getting rid of distance point to own point
in = in(1:mnn, :);
connect = zeros(n, 1);
NN = zeros(mnn, n);

for i = 1:n
    
    % If all connections already attributed, moving on
    con = connect(i);
    if con == mnn
        continue
    end
    
    % Checking for obstacles
    T = zeros(size(in, 1), 1);
    for j = 1:length(T)
        val = in(j, i);
        obstacle = 0;
        for k1 = 1:size(polygons, 1)
            A = polygons{k1};
            ttt = size(A, 1);
            for k2 = 1:2:ttt
                if isLine(pts(i, :), pts(val, :), A(k2, :), A(k2+1, :)) == 1
                    obstacle = 1;
                end
            end
        end
        if obstacle == 0
            T(j) = 1;
        end
    end
    index = in(T == 1, i);
    
    % Check if a point is already listed
    okindex = 1;
    for j = 1:mnn
        if sum(index == NN(j, i)) ~= 0
            if length(index) == 1
                okindex = 0;
                break
            end
            temptemp = find(index == NN(j, i));
            index = index([1:(temptemp-1), (temptemp+1):end]);
        end
    end
    if okindex == 0
        continue
    end
    
    % Check if a point already has its max connection
    okconnect = 1;
    indj = ones(length(index), 1);
    for j = 1:length(index)
        if connect(index(j)) == mnn
            if length(index) == 1
                okconnect = 0;
                break
            end
            indj(j) = 0;
        end
    end
    if okconnect == 0
        continue
    end
    index = index(indj == 1);
    
    % Adding connections
    mnni = length(index);
    if mnni == 0
        continue
    end
    reste = min([mnni, size(NN, 1)-con]);
    NN((con+1):(con+reste), i) = index(1:reste);
    
    % Updating connections
    connect(i) = connect(i) + mnni;
    for j = 1:reste
        valj = index(j);
        if connect(valj) == mnn
            continue
        end
        NN((connect(valj)+1), valj) = i;
        connect(valj) = connect(valj) + 1;
    end
end
        
% Plotting lines
for i = 1:mnn
    for j = 1:n
        if NN(i, j) ~= 0
            value = NN(i, j);
            line([pts(j, 1), pts(value, 1)], [pts(j, 2), pts(value, 2)])
            hold on
            k = find(NN(:, value) == j);
            NN(k, value) = 0;
        end
    end
end



    
    
%% Functions  
function lineornot = isLine(p1, p2, pline1, pline2)
    % Function that will return 1 if there is a line between two points,
    % Returns 0 otherwise.
    
    % Take a point just near pline1 following the vector pline1pline2
    ptest1 = pline1 + 0.00001 * (pline2 - pline1);
    % Take a point just near pline2 following the vector pline2pline1
    ptest2 = pline2 + 0.00001 * (pline1 - pline2);
    % Checking if ptest1 in polynom created by pline1, p1 & p2
    inpoly1 = inpolygon(ptest1(1), ptest1(2), [pline1(1), p1(1), p2(1)], ...
                    [pline1(2), p1(2), p2(2)]);
    % Checking if ptest2 in polynom created by pline2, p1 & p2
    inpoly2 = inpolygon(ptest2(1), ptest2(2), [pline2(1), p1(1), p2(1)], ...
                    [pline2(2), p1(2), p2(2)]);
                
    lineornot = (inpoly1 == 1 & inpoly2 == 1);  
end