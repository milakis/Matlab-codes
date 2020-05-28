%% Code description: 
% Solution PS2b - Exercise 3, Question 6
% Policy function iteration - PFI using using the endogenous grid method
% International Economics and Finance - University of Bonn (SoSe 2020)
% Professor: Keith Kuester

%% Housekeeping 
clear       % clear Workspace
close all   % close all figures
clc         % clear Command Window

%% Start stopwatch timer
tic

%% Parameters
rpar   = 0.02;      % international fixed real rate of interest
betpar = 0.975;     % time discount factor
sigpar = 4;         % coefficient of relative risk aversion

%% Markov transition matrix
p_LL = 0.95;                            % persistence probability for low income
p_MM = 0.95;                            % persistence probability for medium income
p_HH = 0.95;                            % persistence probability for high income
p_MH = 0.05;                            % transition probability from high to medium income
p_ML = 0.05;                            % transition probability from low to medium income
p_LM = 0.025;                           % transition probability from medium to low income
p_HM = 0.025;                           % transition probability from medium to low income
p_LH = 1 - p_HH - p_MH;                 % transition probability from high to low income
p_HL = 1 - p_LL - p_ML;                 % transition probability from low to high income
TransYmat = [p_LL p_ML, p_HL;           % Markov transition probabilities matrix
             p_LM p_MM p_HM;
             p_LH p_MH p_HH];
       
%% Asset grid
bmin = -2;                          % minimum asset level (borrowing constraint)
bmax = 5;                           % maximum asset level (saving constraint)
nasset = 300;                       % number of grid points
stepsize = (bmax-bmin)/(nasset-1);  % distance between grid points (evenly-spaced grid)
bgrid = [bmin:stepsize:bmax]';      % asset grid vector

%% Incomes and assets matrix
y_L = 0.85;                         % income in low state
y_M = 1;                            % income in medium state
y_H = 1.15;                         % income in high state
yvec = [y_L;y_M;y_H];               % vector of income states
ny = length(yvec);                  % number of income states 
Ymat = kron(yvec,ones(1,nasset));   % income matrix
Bmat = kron(ones(ny,1), bgrid');    % asset grid matrix

%% Maximum debt level that can be sustained for sure
dmax = 1/rpar*min(yvec); 
if bmin<-dmax
    error('borrowing limit not tight enough')
end

%% Average income
% Compute unconditional probabilities for states
[eigvectors,eigvalues] = eig(TransYmat','vector');  % eigenvectors and eigenvalues of Markov transition matrix
pos = find(eigvalues==1);                           % position of eigenvalue equal to one
% Note: a Markov matrix always has an eigenvalue 1. 
% All other eigenvalues are in absolute value smaller or equal to 1.
eigvec = eigvectors(:,pos);                         % eigenvector corresponding to to eigenvalue equal to one
probuncondvec = eigvec/sum(eigvec);                 % probabilities for each state
% Compute average income
yavg = probuncondvec'*yvec; 

%% Initial guess and initialization
Cpol = Ymat + rpar*Bmat;      % initial guess for consumption policy function
Cpoln = 0*Cpol;               % initialize consumption matrix
Bpast = 0*Cpol;               % initalize past assets matrix

%% Policy function iteration (PFI) - iteration on the Euler equation
difference = 1.e8;                  % initialize convergence criterion to a big number
convcrit = 1.e-8;                   % convergence criterion for PFI (tolerance level)
counter = 0;                        % iteration counter
maxiter = 1000;                     % maximum number of iterations
confidx = 1:floor(nasset*0.6);      % which parts of the consumption policy function to check for convergence

while ((counter<maxiter) && (difference>convcrit))
    counter  = counter + 1;  %  update counter
    Cpol_old = Cpol; 
    for inc = 1:ny % loop over current income state
        yhere = yvec(inc);
        for asset = 1:nasset % loop over tomorrow's debt level as state
            bhere = bgrid(asset);  
            % consumption in t+1 (using the policy function)
            cplus = Cpol(:,asset);
            % marginal utility of consumption in t+1
            ucplus = cplus.^(-sigpar);
            % expected marginal utility of consumtpion in t+1 as of t
            EMUhere = TransYmat(inc,:)*ucplus;
            % RHS of Euler equation
            RHShere = betpar*(1+rpar)*EMUhere;
            % invert to get current consumption 
            chere = RHShere^(-1/sigpar);
            % invert to get asset at beginning of period
            bpast = (chere+bhere-yhere)/(1+rpar);
            % update consumption and assets matrices
            Cpoln(inc,asset) = chere;
            Bpast(inc,asset) = bpast;
        end
        
        % check if borrowing constraint is binding
        idx_bc = find(bgrid <= Bpast(inc,1));
        % for these points, HH goes to borrowing constraint
        % update consumption function accordingly
        if ~isempty(idx_bc)
            Cpol(inc,idx_bc) = yhere - bgrid(1) + bgrid(idx_bc)*(1+rpar);
        else
            idx_bc = 0;
        end
        % linear interpolation if borrowing constraint is not binding
        cpolvec = interp1(Bpast(inc,:), Cpoln(inc,:), bgrid(max(idx_bc)+1:end),'linear','extrap')';        
        Cpol(inc,max(idx_bc)+1:end) = cpolvec;
        
        % check if savings constraint is binding
        idx_sav = find(bgrid >= Bpast(inc,end));
        % for these points, HH goes to saving constraint
        if ~isempty(idx_sav)
            Cpol(inc,idx_sav) = yhere - bgrid(end) + bgrid(idx_sav)*(1+rpar);
        else
            idx_sav = length(bgrid)+1;
        end
    end
    
    % update consumption policy function
    difference = max(max(abs(Cpol_old(:, confidx)-Cpol(:, confidx))))
    if difference<convcrit
        fprintf('The policy function has converged in %2.0f iterations.\n',counter);
        difference = 0;
    end
    
end

%% Other policy functions
Bpol = -Cpol + (1+rpar)*Bmat + Ymat;    % asset policy function
Savpol = Bpol-Bmat;                     % net savings policy function

%% Plot policy functions
f1 = figure(1);
f1.WindowState = 'maximized';
% Asset policy function
subplot(1,3,1); 
plot(bgrid, Bpol, 'LineWidth', 2); 
legend('y_L','y_M','y_H', 'Location', 'northwest')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('tomorrow''s asset level', 'FontSize', 14)
% Savings policy function
subplot(1,3,2);
plot(bgrid, Savpol, 'LineWidth', 2); 
legend('y_L','y_M','y_H', 'Location', 'northeast')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('today''s savings', 'FontSize', 14)
ylim([-0.4 0.6])
% Consumption policy function
subplot(1,3,3);
plot(bgrid, Cpol, 'LineWidth', 2); 
legend('y_L','y_M','y_H', 'Location', 'northwest')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('consumption', 'FontSize', 14)
ylim([0.4 2])

print('3q6','-dpng','-r300');

%% End stopwatch timer
toc
