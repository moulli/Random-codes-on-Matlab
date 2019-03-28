clear; close all; clc





%% Definition of the parameters:
po = [1, 9]; %initial point
pf = [9, 9]; %final point
large = 10; %window wideness
long = 10; %window height
dmax = 0.3; %maximum distance between nodes
pro = 0.2; %probability to fire towards the final point (goal bias)
imax = 10000; %number of iterations max
sp = 0.0002; %plotting speed
threshobst = 30; %number of time a node blocked due to obstacle

% RRT*
star = 1; %if 1, activate the RRT*
dx = 0.2;


%% Definition of the obstacles
lines = {[3, 3], [3, 10];
         [3, 10], [4, 10];
         [4, 10], [4, 2];
         [4, 2], [2, 2];
         [2, 2], [2, 3];
         [2, 3], [3, 3];
         [5.5, 0], [5.5, 8];
         [5.5, 8], [6, 8];
         [6, 8], [6, 0];
         [6, 0], [5.5, 0]
         [7.5, 3], [7.5, 10];
         [7.5, 10], [8.5, 10];
         [8.5, 10], [8.5, 3];
         [8.5, 3], [7.5, 3]};

     
%% Definition of the window
ax = [0 large 0 long];
msize = 5;
scatter(po(1), po(2), msize+5)
scatter(pf(1), pf(2), msize+5)
title("Rapidly-exploring Random Tree", "Interpreter", "latex")
axis(ax)
hold on
for i = 1:size(lines, 1)
    p1 = lines{i, 1};
    p2 = lines{i, 2};
    line([p1(1), p2(1)], [p1(2), p2(2)], "Color", "black", "Linewidth", 2)
    hold on
end


%% Obtaining the path back
% To obtain the path from the final to the first point, we are going to
% keep record of all the paths, using a vector linking each point to its
% generation point
path = [1]; %first point arbitrarly comes from first point


%% Tree algorithm
P = [po]; %keeping record of points
Obs = [0]; %number of time a node blocked
Pdist = [po]; %number of nodes kept for the distances
success = 0; %did we reach the objective
for i = 1:imax
    % Time sleep
    pause(sp);

    % Deactivating if node is blocked by obstacle
    for k = 1:length(Obs)
        if Obs(k) >= threshobst;
            Pdist(k) = NaN;
        end
    end    
    
    % Random point with probability
    pri = rand(1);
    if pri <= pro
        p_i = pf;
    else
        p_i = [large, long] .* rand(1, 2);
    end
    % Computing distance from all existing nodes
    dist = sqrt(sum((Pdist - p_i).^2, 2));
    % RRT or RRT*
    [dmin, in] = sort(dist);
    Td = dmin <= (dmax + dx); %nodes within a set distance
    if star == 0 | Td == 0
        % Closest node
        dmin = dmin(1);
        in = in(1);
    else
        in = in(Td);
        dmin = dmin(Td);
        pathin = zeros(length(in), 1);
        for j = 1:length(in)
            fpathin = [in(j)];
            while fpathin(end) ~= 1
                fpathin = [fpathin; path(fpathin(end))];
            end
            pathin(j) = length(fpathin);
        end
        [valmin, indmin] = min(pathin); %getting fastest path
        in = in(indmin);
        dmin = dmin(indmin);
        % Problem with last distance if p_i = pf. Let us settle this:
        if p_i == pf & dmin <= dmax + dx
            dx = dx - 0.1;
        end
    end    
    pmin = P(in, :);
    % Changing distance if bigger than threshold
    if dmin > dmax
        p_i = ((dmax / dmin) * (p_i - pmin)) + pmin;
    end
    % Checking p_i does not already exist
    exis = P == p_i;
    if sum(exis) ~= 0
        continue
    end
    
    % If obstacle, then we ignore and keep going
    obstacle = 0; %setting back obstacle to 0
    for j = 1:size(lines, 1)
        p1 = lines{j, 1};
        p2 = lines{j, 2};
        if isLine(pmin, p_i, p1, p2) == 1   
            obstacle = 1;
        end
    end
    if obstacle == 1
        % We are going to try something
        % In certain situations with a goal bias, certain nodes are worthless
        % because they always produce a nonsense due to obstacles
        % We are going to set a threshold for each node, then deactivate it
        Obs(in) = Obs(in) + 1;
        continue %start another loop to avoid obstacle
    end
  
    % Adding point to nodes
    P = [P; p_i];
    Obs = [Obs; 0]; %adding the obstacle parameter
    Pdist = [Pdist; p_i]; %and node to nodes kept for distances
    % Keeping track of the connexions
    path = [path; in]; %new point comes from this point
    % Plotting
    scatter(p_i(1), p_i(2), msize)
    hold on
    line([pmin(1), p_i(1)], [pmin(2), p_i(2)])
    hold on
    % Stopping algorithm if really close to final point
    if sqrt(sum((pf - p_i).^2, 2)) <= 0.001
        success = 1;
        break
    end
end


%% Plotting final path
pause(0.5)
if success == 1
    final_path = length(path);
else
    [Ptop, itop] = min(sqrt(sum((P - pf).^2, 2)));
    final_path = itop;
end
while final_path(end) ~= 1
    final_path = [final_path; path(final_path(end))];
end
Matpath = P(final_path, :);
plot(Matpath(:, 1), Matpath(:, 2), "r", 'LineWidth', 2)
hold on


%% Smoothing path
Matpath = flip(Matpath);
psmooth = [po];
dep = 1;
for i = 1:(length(Matpath)-1)
    for j = 1:size(lines, 1)
        p1 = lines{j, 1};
        p2 = lines{j, 2};
        if isLine(Matpath(dep, :), Matpath(i+1, :), p1, p2) == 1   
            psmooth = [psmooth; Matpath(i, :)];
            dep = i;
            break
        end
    end
end
psmooth = [psmooth; pf];
plot(psmooth(:, 1), psmooth(:, 2), "-.g", 'LineWidth', 2)



    
    


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