clear; close all; clc



% The field is a 5x5 square.
% The agent is characterized by 4 parameters (input vector):
% - its position (it can move or stay in place at each iteration),
% - its hunger (ranges from 0 to 20, dies at 20/ food gives 5 back/ if 
% agent eats while hunger at 0, it dies/ +1 per turn),
% - its tireness (ranges from 0 to 20, dies at 20/ can rest instead of
% moving during an iteration, and lose 2 points of tireness/ +1 per turn),
% - its age (+1 per turn, starts at 0, dies at 100).
% The agent is curious, so it tries to go to the positions it has not been
% in a while (it gets a reward each time it does).
% If agent is on a border and tries to go through, it does not move, but
% is more hungry and less awake.
% The output vector is a number between 1 (up) and 6 (sleep), with:
% 2 (right), 3 (down), 4 (left), and 5 (stay still).
% The agent has a 90% chance of chosing the best option, and 10% chance
% of chosing a random option (exploration).

% Parameters:
epoch = 1000;
dtplot = 0.001;
hmax = 20;
tmax = 20;
food = 19;
alpha = 0.1; % learning rate
gamma = 0.9; % discount factor
rfood = 1; % food reward
rsleep = 1; % sleep reward
rdiscover = 1; % discovering reward
rwall = -1; % hitting a wall reward
rstarve = -1; % hunger gets to 20 reward
rtire = -1; % tireness gets to 20 reward
% Ponderate all these rewards with state of hunger, tireness and discovery!

% Initialization:
input = [1, 0, 0, 0]';
Q = rand(6, length(input));
Qinit = Q;
visited = zeros(25, 1);
visited(1) = 1;
lifetime = [];

% Algorithm:
space = zeros(5);
space(4, 4) = 0.5;
space(input(1)) = 1;
figure
image(space, 'CDataMapping', 'scaled')
colorbar
for i = 1:epoch
    action = Q * input;
    [~, iaction] = max(action);
    % Action or random movement:
    rmov = rand;
    if rmov > 0.9
        iaction = [1:iaction-1, iaction+1:6];
        iaction = iaction(randperm(5, 1));
    end
    % New position:
    ninput = input;
    if iaction == 1 && mod(input(1), 5) ~= 1
        ninput(1) = input(1)-1;
    elseif iaction == 2 && sum(input(1) == [21, 22, 23, 24, 25]) == 0
        ninput(1) = input(1)+5;
    elseif iaction == 3 && mod(input(1), 5) ~= 0
        ninput(1) = input(1)+1;
    elseif iaction == 4 && sum(input(1) == [1, 2, 3, 4, 5]) == 0
        ninput(1) = input(1)-5;
    end        
    % Definition of reward:
    if ninput(1) == input(1) || input(2) == 20 || (input(2) == 0 && input(1) == 19) || input(3) == 20 || (input(3) == 0 && iaction == 6)
        reward = -1;
    elseif visited(ninput(1)) == min(visited) || (ninput(1) == 19 && ninput(2) > 10) || (iaction == 6 && ninput(3) > 10)
        reward = 1;
    else
        reward = 0;
    end
    % Definition of discount:
    if (mod(ninput(1), 5) ~= 1 && visited(ninput(1)-1) == min(visited)) || (sum(ninput(1) == [21, 22, 23, 24, 25]) == 0 && visited(ninput(1)+5) == min(visited)) || ...
            (mod(ninput(1), 5) ~= 0 && visited(ninput(1)+1) == min(visited)) || (sum(ninput(1) == [1, 2, 3, 4, 5]) == 0 && visited(ninput(1)-5) == min(visited)) || ...
            (sum(ninput(1)+[-1, 5, 1, -5] == 19) == 1 && ninput(2)+1 > 10) || (ninput(3)+1 >10)
        discount = 1;
    else
        discount = 0;
    end
    % Updating Q:
    Qline = Q(iaction, :);
    Qline = (1-alpha) * Qline + alpha * (reward + gamma*discount);
    Q(iaction, :) = Qline;
    % Plotting:
    pause(dtplot)
    space = zeros(5);
    space(4, 4) = 0.5;
    if iaction == 6
        space(ninput(1)) == 0.6;
    else
        space(ninput(1)) = 1;
    end
    image(space, 'CDataMapping', 'scaled')
    colorbar
    % Updating:
    input = ninput;
    input(2) = input(2)+1;
    input(3) = input(3)+1;
    input(4) = input(4)+1;
    visited(input(1)) = visited(input(1))+1;
    if input(1) == 19
        input(2) = max(0, input(2)-5);
    end
    if iaction == 6
        input(3) = max(0, input(3)-2);
    end
    if input(2) == 20 || input(3) == 20
        lifetime = [lifetime, input(4)];
        input = [1, 0, 0, 0]';
        visited = zeros(25, 1);
        visited(1) = 1;
    end
end

figure
subplot(1, 3, 1)
image(Qinit, 'CDataMapping', 'scaled')
colorbar
subplot(1, 3, 2)
image(Q, 'CDataMapping', 'scaled')
colorbar
subplot(1, 3, 3)
plot(lifetime)

    




