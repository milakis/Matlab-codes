%% Solve the Cake eating problem

clc
clear
close all

%% Define Numerical Parameters
mpar.nk   = 10;  % number of points on the capital grid
mpar.nz   = 2;    % number of points on the productivity grid
mpar.mink = 0.1;  % lowest point on the capital grid
mpar.maxk = 0.4;    % highest point on the capital grid
mpar.crit = 1e-6; % Precision up to which to solve the value function

%% Define Economic parameters
par.beta  = 0.95; % Discount factor
par.alpha = 0.5;    % Curvature of production function
par.gamma = 1;    % Coefficient of relative risk aversion
par.delta = 1;    % Depreciation
prob.z     = [0.875, 0.125; 0.125, 0.875];% Transition probabilities for productivity

%% Produce grids

grid.k = [mpar.mink:(mpar.maxk-mpar.mink)/(mpar.nk-1):mpar.maxk];
grid.z = [0.9,1.1];
%% Display Model
TablePar={'Discount Factor:', par.beta; 'Returns to Scale', par.alpha; ...
    'Relative Risk Aversion', par.gamma; 'Depreciation', par.delta};

TablePar
%% Define utility functions

if par.gamma ==1
    util  = @(c)log(c); %utility function
    mutil = @(c)1/c; %marginal utility function
else
    util  = @(c)c^(1-par.gamma)/(1-par.gamma);
    mutil = @(c)-c^(-par.gamma);
end

%% Calculate Consumption and Utility for Capital Choices
[meshes.k,  meshes.kprime, meshes.z]= ndgrid(grid.k,grid.k,grid.z);
Y = meshes.k.^par.alpha.*meshes.z; %Income/Resources
C = Y - (1-par.delta).*meshes.k - meshes.kprime; %Consumption

U      = util(C); %Dimensions k x k' x z
U(C<0) = -Inf; % Disallow negative consumption

%% Value Function Iteration

V     = zeros(mpar.nk,mpar.nz);
dist  = 9999;
count = 1;
while dist(count)>mpar.crit
    count       = count+1;                % count the number of iterations
    EV          = par.beta.*V*prob.z.';   % Calculate expected continuation value
    EVfull      = repmat(reshape(EV, [1 mpar.nk mpar.nz]),[mpar.nk 1 1]); % Copy Value to second dimension
    Vnew        = max(U+EVfull,[],2);   % Update Value Function
    dist(count) = max(abs(Vnew(:)-V(:))); % Calculate distance between old guess and update
    V           = squeeze(Vnew);          % Copy update to value function
end

%% Produce Policy functions
[~,policy]  = max(U+EVfull,[],2);
kprime      = grid.k(squeeze(policy));

%% Plots

figure(1)
semilogy((dist(2:end)))
title('Distance between two updates of V -logscale')

figure(2)
plot(grid.k,kprime)
hold on
plot(grid.k,grid.k,'k--')
legend({'Capital policy','45 degree line'})
