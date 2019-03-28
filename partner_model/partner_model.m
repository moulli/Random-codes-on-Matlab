clear; close all; clc



%% Model for number of sexual partners depending on the gender.

% GOAL: see how the number of partners increases depending on parameters.
%
% HYPOTHESE:
% - There are a certain number of populations
% - In each population you can chose the number of males and females
% - It's a population model, it's not at the individual scale
%
% ROUTINE FOR EVERY WEEK:
% - Every week, a chosen proportion of each gender of each population (EGEP) goes out
% - There is a certain variance for this, for EGEP -> gaussian probability
% - The number of persons someone is going to approach is modeled by an exponential law
% - The parameter for this law is different for EGEP
% - Now every one going out has a number of people they are going to approach
% - Everyone has a conditional probability of accepting to ken with a
% random person from the other sex, given the other person want to ken them
% - There is a variance for this parameter as well -> Gaussian
% - For each person approached, we compute the probability they are going
% to accept to ken, and is this is higher than a random number, they
% actually ken, and we increase their partner count
%
% HOW MEMORY WORK:
% Whenever you go out, there is a chance you meet someone you already
% kened. If this happen while you only aproached one person, then it
% cancels this approach: this means either you kened the person you already
% kened, either you kened no one. If you have several approaches, then you
% can approach other persons. Everytime you meet someone you already kened,
% it deletes one of your approches, though
%
% PLOTTING:
% - We average the number of partners for EGEP, and plot them in the end
% - The goal is to have the most accurate parameters for the model
%
% Couple, voyages, décès, ajouter une option qui fait qu'on est
% plus curieux au début, puis on fait moins de nouvelles expériences, mais
% en même temps on ken plus facilement




%% Parameters for the model, and different matrices:

nweeks = 52;


% Population matrix:
% (lines: male or female, column: population)
Npop = 50 * ones(2, 2);
% Npop = [200, 200;
%         200, 200];
    
    
% Number of people of each gender each population to go out:
% (lines: male or female, column: population)
Nout = round(Npop./2);
% If matrix does not correspond to number of population, return error:
if size(Nout, 1) ~= size(Npop, 1) || size(Nout, 2) ~= size(Npop, 2)
    error("Nout matrix not suited to populations")
end

% Variance for number of people of each gender each population to go out:
% (lines: male or female, column: population)
Noutv = [round(Npop(1, :)./20);
         round(Npop(2, :)./10)];
% If matrix does not correspond to number of population, return error:
if size(Noutv, 1) ~= size(Npop, 1) || size(Noutv, 2) ~= size(Npop, 2)
    error("Noutv matrix not suited to populations")
end

    
% Conditional probability to ken for each gender of each population:
% (If other wants to ken, what is the probability that I ken)
% (lines: male or female, column: population)
Nken = [0.6, 0.55;
        0.05, 0.03];
% If matrix does not correspond to number of population, return error:
if size(Nken, 1) ~= size(Npop, 1) || size(Nken, 2) ~= size(Npop, 2)
    error("Nken matrix not suited to populations")
end

% Variance to ken for each gender of each population:
% (lines: male or female, column: population)
Nkenv = [0.05, 0.05;
         0.3, 0.2];
% If matrix does not correspond to number of population, return error:
if size(Nkenv, 1) ~= size(Npop, 1) || size(Nkenv, 2) ~= size(Npop, 2)
    error("Nkenv matrix not suited to populations")
end


% Number of people approached each night out for each gender of each population:
% (How many people am I going to approach in one night: exponential law)
% (lines: male or female, column: population)
Napp = [0.8, 0.6;
        0.2, 0.1];
% If matrix does not correspond to number of population, return error:
if size(Napp, 1) ~= size(Npop, 1) || size(Napp, 2) ~= size(Npop, 2)
    error("Napp matrix not suited to populations")
end


% Age of first time:
% (lines: male or female, column: population)
Nfirst = [16.5*52*ones(1, 2);
          17.5*52*ones(1, 2)];
% If matrix does not correspond to number of population, return error:
if size(Nfirst, 1) ~= size(Npop, 1) || size(Nfirst, 2) ~= size(Npop, 2)
    error("Nfirst matrix not suited to populations")
end

% Variance for age of first time:
% (lines: male or female, column: population)
Nfirstv = [3*ones(1, 2);
           3*ones(1, 2)];
% If matrix does not correspond to number of population, return error:
if size(Nfirstv, 1) ~= size(Npop, 1) || size(Nfirstv, 2) ~= size(Npop, 2)
    error("Nfirstv matrix not suited to populations")
end


% Age of death:
% (lines: male or female, column: population)
Ndeath = [70*52*ones(1, 2);
          80*52*ones(1, 2)];
% If matrix does not correspond to number of population, return error:
if size(Ndeath, 1) ~= size(Npop, 1) || size(Ndeath, 2) ~= size(Npop, 2)
    error("Ndeath matrix not suited to populations")
end

% Variance for age of first time:
% (lines: male or female, column: population)
Ndeathv = [10*ones(1, 2);
           10*ones(1, 2)];
% If matrix does not correspond to number of population, return error:
if size(Ndeathv, 1) ~= size(Npop, 1) || size(Ndeathv, 2) ~= size(Npop, 2)
    error("Ndeathv matrix not suited to populations")
end


%% Creation of the global recap matrix:

V = zeros(8, sum(Npop(:)));
% 1 for population
% 2 for gender
% 3 for actual age
% 4 for age of first time
% 5 for age of death
% 6 for population in which holidays
% 7 for number of weeks left of holidays
% 8 for number of partners
for i = 1:size(Npop, 2)
    V(1, sum(sum(Npop(:, 1:i-1)))+1:sum(sum(Npop(:, 1:i)))) = i; % population
    V(2, sum(sum(Npop(:, 1:i-1)))+1:sum(sum(Npop(:, 1:i-1)))+Npop(1, i)) = 1; % males
    V(2, sum(sum(Npop(:, 1:i-1)))+Npop(1, i)+1:sum(sum(Npop(:, 1:i)))) = 2; % females
end




%% Creation of the list of partners cell:
% Everytime someone kens, other person saved, so that it does not count as
% a new partner again

LP = cell(sum(Npop(:)), 1);




%% Algorithm, for the number of weeks selected:


% Mean number of partners for each gender each population:
NP = zeros(size(Npop, 1)*size(Npop, 2), nweeks + 1);

% Initialization of people:
for j = 1:size(Npop, 2)
    popAnte = sum(sum(Npop(:, 1:j-1)));
    firstH = Nfirst(1, j) + Nfirstv(1, j)*randn(1, Npop(1, j));
    firstF = Nfirst(2, j) + Nfirstv(2, j)*randn(1, Npop(2, j));
    V(4, popAnte+1:sum(sum(Npop(:, 1:j)))) = [firstH, firstF];
    deathH = Ndeath(1, j) + Ndeathv(1, j)*randn(1, Npop(1, j));
    deathF = Ndeath(2, j) + Ndeathv(2, j)*randn(1, Npop(2, j));
    V(5, popAnte+1:sum(sum(Npop(:, 1:j)))) = [deathH, deathF];
end

for i = 1:nweeks
    
    % Updating age:
    V(3, :) = V(3, :) + 1;
    
    for j = 1:size(Npop, 2)
        
        % Taking part of V that interests us:
        Vj = V(:, sum(sum(Npop(:, 1:j-1)))+1:sum(sum(Npop(:, 1:j))));
        
        % Updating age state (should be outside loop, but timestep is only a week, so it's fine)
        
       
        % Who goes out:
        kout = round(Nout(:, j) + randn(2, 1) .* Noutv(:, j));
        goH = randperm(Npop(1, j));
        goH = goH(1:kout(1)); % list of men going out
        goF = randperm(Npop(2, j));
        goF = goF(1:kout(2)); % list of women going out
        
        
        % Number of people each is going to meet:
        npH = round(exprnd(Napp(1, j), 1, length(goH)));
        npF = round(exprnd(Napp(2, j), 1, length(goF)));
        
        
        % Who is going to ken who:
        
        popAnte = sum(sum(Npop(:, 1:j-1))); % to simplify later
        
        % Women first:
        for k = 1:length(goF)
            while npF(k) > 0 && isempty(goH) == 0
                u_ken = rand < Nken(1, j) + randn * Nkenv(1, j);
                if u_ken == 1 && isempty(find(LP{popAnte+Npop(1, j)+goF(k)} == popAnte+goH(1))) == 1
                    Vj(8, goH(1)) = Vj(8, goH(1)) + 1; % updating number of partners
                    Vj(8, Npop(1, j)+goF(k)) = Vj(8, Npop(1, j)+goF(k)) + 1;
                    LP{popAnte+Npop(1, j)+goF(k)} = [LP{popAnte+Npop(1, j)+goF(k)}; popAnte+goH(1)]; % updating list of partners
                    LP{popAnte+goH(1)} = [LP{popAnte+goH(1)}; popAnte+Npop(1, j)+goF(k)];
                    goH = goH(2:end); % deleting male who kenned
                    npH = npH(2:end); 
                    npF(k) = 0; % to exit loop
                end    
                npF(k) = npF(k) - 1;
            end
        end
        
        % Then men: 
        restF = goF(npF >= 0); % deleting women who already kenned
        for k = 1:length(goH)
            while npH(k) > 0 && isempty(goF) == 0
                u_ken = rand < Nken(2, j) + randn * Nkenv(2, j);
                if u_ken == 1 && isempty(find(LP{popAnte+goH(k)} == popAnte+Npop(1, j)+goF(1))) == 1
                    Vj(8, goH(k)) = Vj(8, goH(k)) + 1; % updating number of partners
                    Vj(8, Npop(1, j)+goF(1)) = Vj(8, Npop(1, j)+goF(1)) + 1;
                    LP{popAnte+Npop(1, j)+goF(1)} = [LP{popAnte+Npop(1, j)+goF(1)}; popAnte+goH(k)]; % updating list of partners
                    LP{popAnte+goH(k)} = [LP{popAnte+goH(k)}; popAnte+Npop(1, j)+goF(1)];
                    goF = goF(2:end); % deleting male who kenned
                    npF = npF(2:end); 
                    npH(k) = 0; % to exit loop
                end  
                npH(k) = npH(k) - 1;
            end
        end
            
        
        % Computing mean number of partners:
        NP((j-1)*2 + 1, i+1) = mean(Vj(8, 1:Npop(1, j)));
        NP(2*j, i+1) = mean(Vj(8, Npop(1, j)+1:Npop(1, j)+Npop(2, j)));
        
        % Injecting Vj in V:
        V(:, sum(sum(Npop(:, 1:j-1)))+1:sum(sum(Npop(:, 1:j)))) = Vj;
    end
end



%% Plotting the mean number of partners:

Leg = [];
figure

% Plotting for males:
subplot(1, 2, 1)
leg = cell(size(Npop, 2), 1);
for j = 1:size(Npop, 2)
    plot(NP((j-1)*2 + 1, :))
    leg{j} = num2str(j);
    hold on
end
title("Number of partner for males of different populations", "Interpreter", "latex")
xlabel("Number of weeks", "Interpreter", "latex")
ylabel("Mean number of partners", "Interpreter", "latex")
legend(leg, "Location", "best")
grid on

% Plotting for females:
subplot(1, 2, 2)
leg = cell(size(Npop, 2), 1);
for j = 1:size(Npop, 2)
    plot(NP(2*j, :))
    leg{j} = num2str(j);
    hold on
end
title("Number of partner for females of different populations", "Interpreter", "latex")
xlabel("Number of weeks", "Interpreter", "latex")
ylabel("Mean number of partners", "Interpreter", "latex")
legend(leg, "Location", "best")
grid on

    

        









