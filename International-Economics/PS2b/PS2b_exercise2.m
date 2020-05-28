%% Code description: 
% Solution PS2b - Exercise 2 (Value function iteration - VFI)
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
p_LH = 1 - p_LL;                        % transition probability from low to high income
p_HH = 0.95;                            % persistence probability for high income
p_HL = 1 - p_HH;                        % transition probability from high to low income
TransYmat = [p_LL p_LH;p_HL p_HH];      % Markov transition probabilities matrix

%% Asset grid (discretize asset space)
bmin = -2;                              % minimum asset level (borrowing constraint)
bmax = 20;                              % maximum asset level (upper bound on the grid)
nasset = 2000;                          % number of grid points
stepsize = (bmax-bmin)/(nasset-1);      % distance between grid points (evenly-spaced grid)
bvec = [bmin:stepsize:bmax]';           % asset grid vector

%% Incomes matrix
y_L = 0.85;                             % income in low state
y_H = 1.15;                             % income in high state
yvec = [y_L;y_H];                       % vector of income states
ny = length(yvec);                      % number of income states 
Ymat = kron(ones(nasset,1), yvec');     % income matrix

%% Initial guesses
Bpol = kron([1;1],[1:nasset]); % initial guess for asset policy function
Vval = zeros(ny,nasset);       % initial guess for value function

%% Value function iteration (VFI)
difference = 1.e8;             % initialize convergence criterion to a big number
Vn = zeros(ny,nasset);         % initialize value function update
Cpol = kron([1;1],[1:nasset]); % initialize consumption policy function
convcrit = 1.e-8;              % convergence criterion for VFI (tolerance level)
counter = 0;                   % iteration counter
maxiter = 1000;                % maximum number of iterations

while ((counter<maxiter) && (difference>convcrit))
    counter = counter + 1;
    for inc = 1:ny % loop over current income state
        yhere = yvec(inc);
        for asset = 1:nasset % loop over current asset
            bhere = bvec(asset);
            % try debt on the grid
            chere = yhere + bhere*(1+rpar) - bvec; % budget constraint
            % compute expected value in the Bellman equation
            EVhere = TransYmat(inc,:)*Vval(:,:);
            % RHS of Bellman equation
            % add (1-betpar) because of numerical issues when betpar close to 1
            RHShere =(1-betpar)*chere.^(1-sigpar)./(1-sigpar) + betpar*EVhere'; 
            % punish for negative current consumption: ensure that solution respects feasibility constraints
            idx = find(chere<=0);
            RHShere(idx) = -inf;
            % find max of RHS
            [maxval,idx] = max(RHShere);
            % update for value function
            Vn(inc,asset) = maxval;
            % update policy functions
            Bpol(inc,asset) = idx; 
            Cpol(inc,asset) = chere(idx);
            if size(Bpol,1)>2
                keyboard
            end
        end
    end
    % check for convergence
    difference = max(max(abs(Vn-Vval)))
    if difference<convcrit
        fprintf('The value function has converged in %2.0f iterations.\n',counter);
        difference = 0;
    end
    % if no convergence, update value function
    Vval = Vn; 
end

%% Plot policy functions
f = figure;
f.WindowState = 'maximized';

% Asset policy function
subplot(1,3,1)
plot(bvec, bvec(Bpol), 'LineWidth', 2)
legend('y_L','y_H', 'Location', 'northwest')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('tomorrow''s asset level', 'FontSize', 14)
% Savings policy function
subplot(1,3,2)
Savpol = bvec(Bpol)-[bvec'; bvec'];
plot(bvec(1:end-1), Savpol(:,1:end-1), 'LineWidth', 2)
legend('y_L','y_H', 'Location', 'northeast')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('today''s savings', 'FontSize', 14)
ylim([-0.4 0.6])
% Consumption policy function
subplot(1,3,3)
plot(bvec, Cpol, 'LineWidth', 2)
legend('y_L','y_H', 'Location', 'northwest')
xlabel('today''s asset level', 'FontSize', 14)
ylabel('consumption', 'FontSize', 14)
ylim([0.4 2])

print('2q4','-dpng','-r300');

%% End stopwatch timer
toc

