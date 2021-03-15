// --------------------------------------------------------------------
//  C A L I B R A T E  O M E G A
// --------------------------------------------------------------------
@#if incompleteMarkets
    Omega = 0.6000 ;
@#else
    Omega = 0;
@#endif

// --------------------------------------------------------------------
//  H A R D W I R E D  C A L I B R A T I O N 
//   *** (think twice before modifying) ***
// --------------------------------------------------------------------
delta    = 0.0150;
theta    = 6.0000;
kv_y     = 0.0100;
ceu_ce   = 1-(1-0.79)/0.6;
reprat   = 0.5000;
upsilon  = 0.5;
chi      = 0.5000;
lambdaHP = 1600.00;
UPSILON  = 0.5;
phi      = 0.001;
// --------------------------------------------------------------------
//  S A M P L E  M E A N S
// --------------------------------------------------------------------
load estimationData;
sd      = find(date==1983.00);
fd      = find(date==2007.50);
[cstR0,cstY0,cstI0,cstP0,csts0,cstf0,cstW0,cstdf0,cstds0,csts600] = computeMeans(sd,fd);
cstY   = cstY0; 
cstR   = cstR0; 
cstf   = cstf0; 
csts   = csts0; 
cstP   = cstP0; 
cstW   = cstW0; 
cstdf  = cstdf0; 
cstds  = cstds0; 
csts60 = csts600;

sig    = 1.0000;
ni     = 2.0000;
nu     = 0.5000;
hF     = 0.7500;
alpha  = 0.7500;
gammap = 0.5000;
gammaw = 0.5000;
psin   = -0.5000;
rhoTR  = 0.7500;
api    = 1.8000;
ay     = 0.1250;
