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
life = 100;
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
sizeNN = [50]; % size of the neural net
nite = 2000;
batch = 100;
step = 0.01;
lambda = 0;

% Initialization:
input = [1, 0, 0, 0]';
sizeNN = [length(input), sizeNN, 6];
Q = cell(length(sizeNN)-1, 1);
for i = 1:length(sizeNN)-1
    Q{i} = 0.012 * (rand(sizeNN(i+1), sizeNN(i)+1) - 0.5);
end
Qinit = Q;
visited = zeros(25, 1);
visited(1) = 1;
lifetime = [];
decisiontot = [];

% Algorithm:
for k1 = 1:nite
    Xcell = cell(batch, 2);
    for k2 = 1:batch
        Xcell{k2, 1} = [];
        input = [1, 0, 0, 0]';
        visited = zeros(25, 1);
        visited(1) = 1;
        for i = 1:life
            Xcell{k2, 1} = [Xcell{k2, 1}; input'];
            iaction = fForwardProp(Q, input);
%             decisiontot = [decisiontot; iaction];
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
            % Updating:
            input = ninput;
            input(2) = input(2)+1;
            input(3) = input(3)+1;
            input(4) = input(4)+1;
            visited(input(1)) = visited(input(1))+1;
            if input(2) == hmax || input(3) == tmax
                lifetime = [lifetime, input(4)];
                Xcell{k2, 2} = -1;
                break
            elseif input(2) <= 4 && input(1) == 19
                lifetime = [lifetime, input(4)];
                Xcell{k2, 2} = -0.5;
                break
            end
            if input(4) == life && nnz(visited) == length(visited)
                lifetime = [lifetime, input(4)];
                Xcell{k2, 2} = 1;
                break
            elseif input(4) == life && nnz(visited) < length(visited)
                lifetime = [lifetime, input(4)];
                Xcell{k2, 2} = 0.5;
                break
            end
            if input(1) == 19
                input(2) = max(0, input(2)-5);
            end
            if iaction == 6
                input(3) = max(0, input(3)-2);
            end
        end
    end
    Q = fReinforcement(Q, Xcell, step, lambda);
    if mod(k1, floor(nite/10)) == 0
        k1
    end
end

figure
plot(lifetime)
hold on
plot(decisiontot)

    

    




