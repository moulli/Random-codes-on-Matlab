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

nweeks = 42*52;


% We are going to treat the population matrices as classes:

Npop(1).m = [5; % number of males and females in population
             2.5; % number that goes out every week
             0.25; % variance for people going out
             0.55; % conditional probability to ken
             0.05; % variance for conditional probability to ken
             0.6; % exponential law parameter for the approaches
             16.5*52; % mean age of first intercourse
             3*52; % variance for age of first intercourse
             70*52; % mean age of death
             10*52]; % variance for mean age of death
Npop(1).f = [5; 2.5; 0.5; 0.03; 0.2; 0.1; 17.5*52; 3*52; 80*52; 10*52];

Npop(2).m = [5; 2.5; 0.25; 0.6; 0.05; 0.8; 16.5*52; 3*52; 70*52; 10*52];
Npop(2).f = [5; 2.5; 0.5; 0.05; 0.3; 0.2; 17.5*52; 3*52; 80*52; 10*52];

npop = length(Npop); % PUT HERE NUMBER OF POPULATIONS!




%% Creation of the global recap cell:

V = cell(npop, 2);
id = 1; % give each person a special id to recognize them easily

for j = 1:npop
    
    % Initializing males:
    nm = Npop(j).m(1);
    V{j, 1} = zeros(8, nm);
    V{j, 1}(1, :) = j;
    V{j, 1}(2, :) = id:(id+nm-1);
    id = id + nm;
    V{j, 1}(4, :) = Npop(j).m(7) + Npop(j).m(8) * randn(1, nm);
    V{j, 1}(5, :) = Npop(j).m(9) + Npop(j).m(10) * randn(1, nm);
    
    % Initializing females:
    nf = Npop(j).f(1);
    V{j, 2} = zeros(8, nf);
    V{j, 2}(1, :) = j;
    V{j, 2}(2, :) = id:(id+nf-1);
    id = id + nf;
    if j == npop
        id = id -1;
    end
    V{j, 2}(4, :) = Npop(j).f(7) + Npop(j).f(8) * randn(1, nf);
    V{j, 2}(5, :) = Npop(j).f(9) + Npop(j).f(10) * randn(1, nf);
    
end




%% Creation of the list of persons cell:

LP = zeros(id); % matrix to know who kenned who: males rows, females columns




%% Algorithm, for the number of weeks selected:


% Mean number of partners for each gender each population:
NP = zeros(2*npop, nweeks + 1);

% Algorithm for each week:
for i = 1:nweeks
    
    for j = 1:size(Npop, 2)
        
        % Updating age (timestep only a week so it works fine for death):
        % Males:
        V{j, 1}(3, :) = V{j, 1}(3, :) + 1;
        born = find(V{j, 1}(3, :) >= V{j, 1}(5, :));
        n_born = length(born);
        V{j, 1}(3, born) = 0;
        V{j, 1}(4, born) = Npop(j).m(7) + Npop(j).m(8) * randn(1, n_born);
        V{j, 1}(5, born) = Npop(j).m(9) + Npop(j).m(10) * randn(1, n_born);
        V{j, 1}(8, born) = 0;
        id_born = V{j, 1}(2, born);
        LP(id_born, :) = 0;
        % Females:
        V{j, 2}(3, :) = V{j, 2}(3, :) + 1;
        born = find(V{j, 2}(3, :) >= V{j, 2}(5, :));
        n_born = length(born);
        V{j, 2}(3, born) = 0;
        V{j, 2}(4, born) = Npop(j).f(7) + Npop(j).f(8) * randn(1, n_born);
        V{j, 2}(5, born) = Npop(j).f(9) + Npop(j).f(10) * randn(1, n_born);
        V{j, 2}(8, born) = 0;
        id_born = V{j, 2}(2, born);
        LP(:, id_born) = 0;

       
        % Who goes out:
        % Males:
        koutm = round(Npop(j).m(2) + Npop(j).m(3)*randn);
        permM = randperm(size(V{j, 1}, 2));
        goM = V{j, 1}(2, permM);
        goM = goM(1:koutm);
        % Females:
        koutf = round(Npop(j).f(2) + Npop(j).f(3)*randn);
        permF = randperm(size(V{j, 2}, 2));
        goF = V{j, 2}(2, permF);
        goF = goF(1:koutf);
        
        
        % Number of people each is going to meet:
        % Males:
        npM = round(exprnd(Npop(j).m(6), 1, length(goM)));
        npF = round(exprnd(Npop(j).f(6), 1, length(goF)));
        
        
        % Who is going to ken who:
        % Women first:
        for k = 1:length(goF)
            while npF(k) > 0 && isempty(goM) == 0
                u_ken = rand < Npop(j).f(4) + randn * Npop(j).f(5);
                if u_ken == 1 && LP(goM(1), goF(k))
                    V{j, 1}(8, permM(1)) = V{j, 1}(8, permM(1)) + 1; % updating number of partners
                    V{j, 2}(8, permF(k)) = V{j, 2}(8, permF(k)) + 1;
                    LP(goM(1), goF(k)) = 1; % updating list of partners
                    LP(goF(k), goM(1)) = 1;
                    goM = goM(2:end); % deleting male who kenned
                    npM = npM(2:end); 
                    npF(k) = 0; % to exit loop
                end    
                npF(k) = npF(k) - 1;
            end
        end
        % Then men: 
        restF = goF(npF >= 0); % deleting women who already kenned
        for k = 1:length(goM)
            while npM(k) > 0 && isempty(goF) == 0
                u_ken = rand < Npop(j).m(4) + randn * Npop(j).m(5);
                if u_ken == 1 && LP(goM(k), goF(1)) == 0
                    V{j, 1}(8, permM(k)) = V{j, 1}(8, permM(k)) + 1; % updating number of partners
                    V{j, 2}(8, permF(1)) = V{j, 2}(8, permF(1)) + 1;
                    LP(goM(k), goF(1)) = 1; % updating list of partners 
                    LP(goF(1), goM(k)) = 1
                    goF = goF(2:end); % deleting female who kenned
                    npF = npF(2:end); 
                    npM(k) = 0; % to exit loop
                end  
                npM(k) = npM(k) - 1;
            end
        end
            
        
        % Computing mean number of partners:
        NP((j-1)*2 + 1, i+1) = mean(V{j, 1}(8, :));
        NP(2*j, i+1) = mean(V{j, 2}(8, :));
    end
end



%% Plotting the mean number of partners:

Leg = [];
figure

% Plotting for males:
subplot(1, 2, 1)
leg = cell(npop, 1);
for j = 1:npop
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
leg = cell(npop, 1);
for j = 1:npop
    plot(NP(2*j, :))
    leg{j} = num2str(j);
    hold on
end
title("Number of partner for females of different populations", "Interpreter", "latex")
xlabel("Number of weeks", "Interpreter", "latex")
ylabel("Mean number of partners", "Interpreter", "latex")
legend(leg, "Location", "best")
grid on









