%% Solve the Consumption Savings Model using an Endogenous Grid Method (EGM)

clear
clc
close all
%% 1. Define parameters

% Numerical parameters
mpar.nk   = 100;   % Number of points on the asset grid
mpar.nz   = 2;    % Number of points on the log-productivity grid
mpar.crit = 1e-5; % Numerical precision
mpar.maxk = 10;    % Maximimum assets
mpar.mink = 10/3;    % Minimum Assets (equal to Borrowing Limit)
disp('Numerical parameters')
mpar % Display numerical parameters
% Economic Parameters
par.r     = 3/90;% Real Rate
par.gamma = 1;    % Coeffcient of relative risk aversion
par.beta  = 0.95; % Discount factor
par.b     = mpar.mink; % Borrowing Limit
disp('Economic parameters')
par % Display economic parameters

%% 2. Generate grids, Meshes and Income
gri.k   = exp(linspace(log(mpar.mink),log(mpar.maxk),mpar.nk)); %Define asset grid on log-linearspaced
prob.z  = [3/5, 2/5; 4/90,  86/90];
P = [3/5, 2/5; 4/90,  86/90];
gri.z   = [1/9, 10/9];
% Meshes of capital and productivity
[meshes.k,  meshes.z] = ndgrid(gri.k,gri.z);
Y = (1+par.r).*meshes.k + meshes.z; % Cash at hand (Labor income plus assets cum dividend)

%% 3. Define utility functions

if par.gamma ==1
    util     = @(c)log(c); % Utility
    mutil    =  @(c)1./c;  % Marginal utility
    invmutil =  @(mutil)1./mutil;% inverse marginal utility
else
    util     = @(c)c^(1-par.gamma)/(1-par.gamma); % Utility
    mutil    = @(c)-c^(-par.gamma); % Marginal utility
    invmutil =  @(mutil)1./mutil; % inverse marginal utility
end

%% 4. Value Function Iteration

tic % Start timer
V    = zeros(mpar.nk,mpar.nz); % Initialize Value Function
distVF = 1; % Initialize Distance
iterVF = 1; % Initialize Iteration count
while distVF(iterVF)>mpar.crit % Value Function iteration loop: until distance is smaller than crit.
    % Update Value Function using off-grid search
    [Vnew,kprime] = VFI_update_spline(V,Y,util,par,mpar,gri,prob); % Optimize given cont' value
    dd            = max(abs(Vnew(:)-V(:))); % Calculate distance between old guess and update

    V             = Vnew; % Update Value Function
    iterVF        = iterVF+1; %Count iterations
    distVF(iterVF)= dd;   % Save distance
end
time(1)=toc; % Save Time used for VFI
%% 5. Plot Policy Functions from VFI

figure(1)
plot(gri.k,kprime) % Plot policy functions
hold on
plot(gri.k,gri.k,'k--') % Add 45 line
title('Policy Function from VFI') % Title and legend of the graph
legend({'low productivity','high productivity','45 degree line'})

%% 6. Endogenous Grid method using linear interpolation
% $$\frac{\partial u}{\partial c}\left[C^*(k',z)\right]=(1+r) \beta E_{z}\left\{\frac{\partialu}{\partial
% c}\left[C(k',z')\right]\right\}$$

tic % Reset timer
C     = (meshes.z  + (1+par.r)*meshes.k); %Initial guess for consumption policy: roll over assets
Cold  = C; % Save old policy
distEG  = 1; % Initialize Distance
iterEG  = 1; % Initialize Iteration count
while distEG>mpar.crit
    C      = EGM(Cold,mutil,invmutil,par,mpar,prob,meshes,gri); % Update consumption policy by EGM
    dd     = max(abs(C(:)-Cold(:))); % Calculate Distance

    Cold   = C; % Replace old policy
    iterEG = iterEG+1; %count iterations
    distEG(iterEG) = dd;
end
[C,Kprimestar] = EGM(C,mutil,invmutil,par,mpar,prob,meshes,gri);
time(2)        = toc; %Time to solve using EGM
%% 7. Plot Policy Functions from Collocation and compare to VFI

figure(2) %Plot Policy Functions from Collocation
plot(gri.k,Kprimestar) % Plot Policy function from Collocation
hold on
plot(gri.k,gri.k,'k--') % Add 45 line
title('Policy Function from EGM') % Title and Legend
legend({'low productivity','high productivity','45 degree line'})

figure(3) %Plot Differences in Policies
plot(gri.k,Kprimestar - kprime)
title('Difference in Policy Function')

%% 8. Compare times of algorithms
disp('Time to solve (VFI, EGM)')
disp(time)
disp('Iterations to solve (VFI, EGM)')
disp([iterVF iterEG])

%% FUNCTIONS (need to be copied to extra files or run as a "Live Script")

%% VF Update
function [Vnew,kprime] = VFI_update_spline(V,Y,util,par,mpar,gri,prob)
V=reshape(V,[mpar.nk,mpar.nz]);
kprime = zeros(size(V));
Vnew   = zeros(size(V));
EV     = par.beta* V* prob.z';   % Calculate expected continuation value

for zz=1:mpar.nz
    ev_int= griddedInterpolant({gri.k},EV(:,zz),'spline');
    for kk=1:mpar.nk
        f             = @(k)(-util(Y(kk,zz)-k)-ev_int(k));
        [kp,v]        = fminbnd(f,par.b,Y(kk,zz));
        Vnew(kk,zz)   = -v;
        kprime(kk,zz) = kp;
    end
end
Vnew=Vnew(:);
end

%% Policy update by EGM
function [C,Kprime] = EGM(C,mutil,invmutil,par,mpar,prob,meshes,gri)
    %% This function iterates forward the consumption policies for the consumption Savings
    % model using the EGM method. C (k x z) is the consumption policy guess. MUTIL and INVMUTIL are
    % the marginal utility functions and its inverse. PAR and MPAR are parameter structures.
    % P is the transition probability matrix. MESHES and GRI are meshes and grids for income
    % (z) and assets (k).
    
    
    mu     = mutil(C); % Calculate marginal utility from c'
    emu    = mu*prob.z';     % Calculate expected marginal utility
    Cstar  = invmutil(par.beta*(1+par.r)*emu);     % Calculate cstar(m',z)
    Kstar  = (Cstar + meshes.k - meshes.z)/(1+par.r); % Calculate mstar(m',z)
    Kprime = meshes.k; % initialze Capital Policy

    for z=1:mpar.nz % For massive problems, this can be done in parallel
        % generate savings function k(z,kstar(k',z))=k'
        Savings     = griddedInterpolant(Kstar(:,z),gri.k,'linear');
        Kprime(:,z) = Savings(gri.k);     % Obtain k'(z,k) by interpolation

    end
    BC         = meshes.k<repmat(Kstar(1,:),mpar.nk,1); % Check Borrowing Constraint
    % Replace Savings for HH saving at BC
    Kprime(BC) = par.b; % Households with the BC flag choose borrowing contraint
    % generate consumption function c(z,k^*(z,k'))
    C          = (1+par.r).*meshes.k - Kprime + meshes.z; %Consumption update from budget constraint
    

end
