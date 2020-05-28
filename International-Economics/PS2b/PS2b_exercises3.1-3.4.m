%% Code description: 
% Solution PS2b - Exercise 3, Questions 1-4
% (Policy function iteration - PFI using using the endogenous grid method)
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
p_LH = 1 - p_LL;                        % transition probability from high to low income
p_HH = 0.95;                            % persistence probability for high income
p_HL = 1 - p_HH;                        % transition probability from low to high income
TransYmat = [p_LL p_HL;                 % Markov transition probabilities matrix
             p_LH p_HH];      

%% Asset grid
bmin = -2;                          % minimum asset level (borrowing constraint)
bmax = 20;                          % maximum asset level
nasset = 300;                       % number of grid points
stepsize = (bmax-bmin)/(nasset-1);  % distance between grid points (evenly-spaced grid)
bgrid = [bmin:stepsize:bmax]';      % asset grid vector

%% Incomes and assets matrix
y_L = 0.85;                         % income in low state
y_H = 1.15;                         % income in high state
yvec = [y_L;y_H];                   % vector of income states
ny = length(yvec);                  % number of income states 
Ymat = kron(yvec,ones(1,nasset));   % income matrix
Bmat = kron(ones(ny,1), bgrid');    % asset grid matrix

%% Maximum debt level that can be sustained for sure
dmax = 1/rpar*min(yvec); 
if bmin<-dmax
    error('borrowing limit not tight enough')
end

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
        idx = find(bgrid <= Bpast(inc,1));
        % for these points, HH goes to borrowing constraint
        % update consumption function accordingly
        if ~isempty(idx)
            Cpol(inc,idx) = yhere - bgrid(1) + bgrid(idx)*(1+rpar);
        else
            idx = 0;
        end
        % linear interpolation if constraint is not binding
        cpolvec = interp1(Bpast(inc,:), Cpoln(inc,:), bgrid(max(idx)+1:end),'linear','extrap')';        
        Cpol(inc,max(idx)+1:end) = cpolvec;
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
legend('y_L','y_H', 'Location', 'northwest')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('tomorrow''s asset level', 'FontSize', 14)
% Savings policy function
subplot(1,3,2);
plot(bgrid, Savpol, 'LineWidth', 2); 
legend('y_L','y_H', 'Location', 'northeast')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('today''s savings', 'FontSize', 14)
ylim([-0.4 0.6])
% Consumption policy function
subplot(1,3,3);
plot(bgrid, Cpol, 'LineWidth', 2); 
legend('y_L','y_H', 'Location', 'northwest')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('consumption', 'FontSize', 14)
ylim([0.4 2])

print('3q2_085_115','-dpng','-r300');

%% End stopwatch timer
toc

%% Model simulation
% Number of draws
ninit = 1000;    % number of draws for burn-in
ndraw = 1000;    % number of draws thereafter
% Initialize variables
cvec     = zeros(ninit+ndraw,1);    % consumption
bvec     = 0*cvec;                  % assets for tomorrow
ysvec    = 0*cvec;                  % income state
bvec(1)  = 0;                       % initial asset state: zero assets
ysvec(1) = 1;                       % initial income state: low-income state
% Draw shocks
shockvec = rand(ninit+ndraw,1);
% Simulate forward
for ast=1:ninit+ndraw
    % current states
    ystate = ysvec(ast);
    bstate = bvec(ast);
    % derive consumption policy for the states
    % linear interpolation because a state may not be on the grid
    cvec(ast,1) = interp1(bgrid, Cpol(ystate,:), bstate, 'linear','extrap')';
    % derive asset policy from the budget constraint
    bvec(ast+1,1) = yvec(ystate) + bstate*(1+rpar) - cvec(ast,1);
    % update output state (for next period)
    if ystate == 1
        if shockvec(ast)>TransYmat(ystate,ystate)
            ystate = 2;
        end
    elseif ystate == 2
        if shockvec(ast)>TransYmat(ystate,ystate)
            ystate = 1;
        end
    end
    ysvec(ast+1) = ystate;
end

% Define variables
tbvec = yvec(ysvec(1:end-1))-cvec; % trade balance
incvec = yvec(ysvec);              % incomes associated with the income state
% Compute moments
disp(['The correlation between trade balance and income is ', num2str(corr(tbvec(ninit+1:end),ysvec(ninit+1:end-1)))])
disp(['Average consumption  ', num2str(mean(cvec(ninit+1:end,:)))])
disp(['Standard dev. consumption  ', num2str(std(log(cvec(ninit+1:end,:)))*100)])
disp(['Average income  ', num2str(mean(incvec(ninit+1:end,:)))])
disp(['Standard dev. income  ', num2str(std(log(incvec(ninit+1:end,:)))*100)])
disp(['Average assets  ', num2str(mean(bvec(ninit+1:end,:)))])

% Plot density functions
f2 = figure(2);
f2.WindowState = 'maximized';
% Assets density function
subplot(1,4,1); 
histogram(bvec(ninit+1:end),30, 'Normalization', 'pdf'); 
xlabel('assets end of period', 'FontSize', 14);
% Consumption density function 
subplot(1,4,2);
histogram(cvec(ninit+1:end),30,'Normalization', 'pdf'); 
xlabel('consumption', 'FontSize', 14);
% Trade balance density function
subplot(1,4,3);
histogram(tbvec(ninit+1:end),30,'Normalization', 'pdf'); 
xlabel('trade balance', 'FontSize', 14);
% Trade balance (conditional on low income) density function
subplot(1,4,4);
histogram(tbvec(ninit+1:end),30,'Normalization', 'pdf'); xlabel('trade balance');
tb2vec = tbvec(ninit+1:end);
ys2vec = ysvec(ninit+1:end-1);
tb2vec(find(ys2vec==1))=[];
histogram(tb2vec,30,'Normalization', 'pdf'); 
xlabel('trade balance cond. on low income', 'FontSize', 14);

print('simul','-dpng','-r300');

