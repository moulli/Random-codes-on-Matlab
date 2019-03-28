clear; close all; clc





%% Definition of the parameters:
po = [1, 9]; %initial point
pf = [9, 9]; %final point
large = 10; %window wideness
long = 10; %window height
pro = 0.2; %probability to fire towards the final point (goal bias)
imax = 10000; %number of iterations max
sp = 0.0002; %plotting speed



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
nbpath = [0]; %path size 


%% Tree algorithm
P = [po]; %keeping record of points
success = 0; %did we reach the objective
for i = 1:imax
    % Time sleep
    pause(sp);
    
    % Random point with probability
    pri = rand(1);
    if pri <= pro
        p_i = pf;
    else
        p_i = [large, long] .* rand(1, 2);
    end
    
    % Checking p_i does not already exist
    exis = P == p_i;
    if sum(exis) ~= 0
        continue
    end
    
    % Checking which is the best node to link it to
    [valpath, indpath] = sort(nbpath);
    obstacle = 1;
    indint = 0;
    while obstacle ~= 0
        if indint == length(indpath)
            break
        end
        obstacle = 0;
        indint = indint + 1;
        in = indpath(indint);
        pmin = P(in, :);
        for j = 1:size(lines, 1)
            p1 = lines{j, 1};
            p2 = lines{j, 2};
            if isLine(pmin, p_i, p1, p2) == 1   
                obstacle = 1;
            end
        end 
    end
    if obstacle == 1
        continue
    end
  
    % Adding point to nodes
    P = [P; p_i];
    % Keeping track of the connexions
    path = [path; in]; %new point comes from this point
    nbpath = [nbpath; (nbpath(in)+1)];
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