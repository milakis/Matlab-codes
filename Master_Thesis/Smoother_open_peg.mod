// SW-like model with precautionary savings, labor market search
//
// Please, refer to the technical appendix for further details
//
// CMRR - Latest version: 27/01/2014


// --------------------------------------------------------------------
//  D E F I N E  M O D E L  T Y P E
// --------------------------------------------------------------------
@#define incompleteMarkets = 1
// Takes value 1 if markets are incomplete (Omega = 0.6)
// Takes value 0 if markets are complete (Omega=0)

// --------------------------------------------------------------------
//  D E C L A R A T I O N  B L O C K
// --------------------------------------------------------------------
@#include "Declaration_open_peg.mod"

// --------------------------------------------------------------------
//  C A L I B R A T I O N  B L O C K
// --------------------------------------------------------------------
@#for name in [ "c", "w", "i", "p", "R", "z", "s", "a", "me", "m", "t"]
    rho@{name} = 0;
    sig@{name} = 0;
@#endfor
@#include "Calibration.mod"
@#include "Parameters.m"


// --------------------------------------------------------------------
//  M O D E L  B L O C K
// --------------------------------------------------------------------
@#include "Model_open_peg.mod"

steady;
check;

stoch_simul(irf=100,order=1,nodisplay,nograph,noprint);
// varobs dI_o dC_o sh60_o dP_o R_o dW_o s_o f_o;
