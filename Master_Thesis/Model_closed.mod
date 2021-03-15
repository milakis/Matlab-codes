// --------------------------------------------------------------------
//  M O D E L  B L O C K
// --------------------------------------------------------------------
model;
// --------------------------------------------------------------------
//
// A - Preamble
//
// --------------------------------------------------------------------
// A.1 - Parameters not contingent on the value found for psi
// --------------------------------------------------------------------
#muz     = 0;
#fss     = cstf/100;
#sss     = csts/100;
#Rss     = cstR/100;
#piss    = 0;
#rss     = (1+Rss)/(1+piss)-1;
#mu      = theta/(theta-1);
#betaF   = exp(sig*muz)*(1+piss)/(1+Rss);
#MFss    = betaF/exp(sig*muz);
#rhoss   = sss/(1-fss);
#rhobar  = log(rhoss/(1-rhoss));
#nss     = fss/(sss+fss);
#mss     = rhoss*nss;
#bF_bW   = sss*(ceu_ce)^(-sig) + (1-sss);
#betaW   = betaF/bF_bW;
#vss     = mss^(1/(1-chi))/(1-(1-rhoss)*nss)^(chi/(1-chi));
#lambdass= mss/vss;
#psi     = 1;  
#alpo    = Omega + (psi)*(1-Omega);
#yss     = exp(-muz)*(alpo)*nss;
#kappav  = kv_y*yss/vss;
#qss     = 1/mu;
#wss     = qss - kappav/lambdass*(1-(1-rhoss)*MFss);
#wbar    = wss*(exp(-muz)/(1+piss))^(gammaw/(gammaw-1));
#tau     = reprat*(1-nss)/nss;

// --------------------------------------------------------------------
//
// B - Model Equations
//
// --------------------------------------------------------------------
// B.1 - Common structure
// --------------------------------------------------------------------
@#for name in [ "IM","HM"]

    // Eq. 1: Real interest rate
    1 + r_@{name} = exp(fic)*(1+R_@{name}(-1))/(1+pi_@{name});

    // Eq. 2: Firm owner's SDF
    MF_@{name} = betaF*exp(-sig*(muz+fiz))*lamF_@{name}/lamF_@{name}(-1);

    // Eq. 3: Firm owner's MU of consumption
    lamF_@{name} = (cF_@{name})^(-sig);

    // Eq. 4: Firm owner Euler equation (bonds)
    1 = MF_@{name}(+1)*(1+r_@{name}(+1))*exp(fip(+1))/exp(fip);

    // Eq. 9: Workers SDF
    MW_@{name} = betaW*exp(-sig*(muz+fiz))*( (1-s_@{name})*lame_@{name} + s_@{name}*lameu_@{name} )/lame_@{name}(-1);

    // Eq. 10: Family members' MU of consumption
    lame_@{name} = (ce_@{name})^(-sig);

    // Eq. 11: MU of consumption for people e->u
    lameu_@{name} = (ceu_@{name})^(-sig);

    // Eq. 13: Mass of people e->u
    neu_@{name} = s_@{name}*n_@{name}(-1);

    // Eq. 14: Mass of people u->u
    nuu_@{name} = 1 - n_@{name} - neu_@{name};

    // Eq. 15:
    (1+pi_@{name})^(1-theta)=(1-alpha)*(1+pi_star_@{name})^(1-theta)+alpha;

    // Eq. 16:
    (1+pi_star_@{name})=mu*(1+pi_@{name})*x1_@{name}/x2_@{name};

    // Eq. 17:
    x1_@{name}=lamF_@{name}*mc_@{name}*y_@{name}+alpha*betaF*(1+pi_@{name}(+1))^theta*x1_@{name}(+1);

    // Eq. 18:
    x2_@{name}=lamF_@{name}*y_@{name}+alpha*betaF*(1+pi_@{name}(+1))^(theta-1)*x2_@{name}(+1);

    // Eq. 19:
    Delta_@{name}=(1-alpha)*(1+pi_star_@{name})^(-theta)*(1+pi_@{name})^theta+alpha*(1+pi_@{name})^theta*Delta_@{name}(-1);
    
    // Eq. 20: Employment dynamics
    n_@{name} = (1-rho_@{name})*n_@{name}(-1) + lambda_@{name}*v_@{name};

    // Eq. 21: Matches
    m_@{name} = exp(fim)*(1-(1-rho_@{name})*n_@{name}(-1))^chi*v_@{name}^(1-chi);

    // Eq. 22: Separation
    s_@{name} = rho_@{name}*(1-f_@{name});

    // Eq. 23: Job finding rate
    f_@{name} = m_@{name}/(1-(1-rho_@{name})*n_@{name}(-1));

    // Eq. 24: Vacancy filling rate
    lambda_@{name} = m_@{name}/v_@{name};

     // Eq. 25: FOC on vacancies
    kappav/lambda_@{name} = (alpo)*(q_@{name} - w_@{name}) 
        + exp(muz+fiz(+1))*MF_@{name}(+1)*(1-rho_@{name}(+1))*kappav/lambda_@{name}(+1);

    // Eq. 26: Aggregate production
    y_@{name} = exp(-muz-fiz)*(alpo)*n_@{name}/Delta_@{name};

    // Eq. 27: Marginal costs
    mc_@{name} = q_@{name}*exp(muz);

    // Eq. 28
    1+pi_@{name} = p_@{name}/p_@{name}(-1);

    // Eq. 30: Zero profit of insurance companies
    tau*w_@{name}*n_@{name} = bu_@{name}*(1-n_@{name});

    // Eq. 31: Wage rate
    w_@{name} = (w_@{name}(-1)*exp(-muz-fiz)/(1+pi_@{name}))^gammaw*(wbar*(n_@{name}/nss)^psin)^(1-gammaw)*exp(fiw);

    // Eq. 32: Taylor rule
    log((1+R_@{name})) = -log(betaF) + (1+rhoTR)*pi_@{name};

    // Eq. 33: Exogenous separation
    rho_@{name} = 1/(1+exp(-rhobar-fis));

    // Eq. 34: Aggregate consumption 
    c_@{name} = (1-Omega)*cF_@{name} + Omega*(n_@{name}*ce_@{name}+nuu_@{name}*cuu_@{name}+neu_@{name}*ceu_@{name});

    // Eq. 35
    y_@{name}   = c_@{name} + kappav*v_@{name} ; 

@#endfor

// --------------------------------------------------------------------
// B.3 - Eqs. specific to IM
// --------------------------------------------------------------------
  // Eq. 5: Workers Euler equation (bonds)
    1 = MW_IM(+1)*(1+r_IM(+1))*exp(fip(+1))/exp(fip);

    // Eq. 6: Consumption of people u->u
    cuu_IM = bu_IM +(1+r_IM)*(exp(fia(-1))-1)*exp(-muz-fiz)-(exp(fia)-1) ;

    // Eq. 7: Budget constraint of workers' family head
    ae_IM + ce_IM = (1-tau)*w_IM + (1+r_IM)*Ae_IM/n_IM ;

    // Eq. 8: Consumption of people e->u
    ceu_IM = bu_IM + (1+r_IM)*ae_IM(-1)*exp(-muz-fiz)-(exp(fia)-1);
    
    // Eq. 29: Net foreign assets
    0 = (1-Omega)*aF_IM + Omega*Ae_IM;

    // Eq. 12: Definition of family's wealth
    exp(muz + fiz)*Ae_IM = (1-s_IM)*n_IM(-1)*ae_IM(-1);
// --------------------------------------------------------------------
// B.4 - Eqs. specific to HM
// --------------------------------------------------------------------
  
   // Eq. Spec 1: Workers Euler equation (bonds)
   ae_HM = 0;

   // Eq. Spec 2: Budget constraint of workers' family head
   ce_HM = w_HM*n_HM ;

   // Eq. Spec 3: Consumption of people e->u
   ceu_HM = w_HM*n_HM;

   // Eq. Spec 4: Consumption of people u->u
   cuu_HM= w_HM*n_HM;
   
    // Eq. 29: Net foreign assets
    aF_HM = 0;

    // Eq. 12: Definition of family's wealth
    Ae_HM = 0;


// --------------------------------------------------------------------
//
// C - Shocks (common to both versions)
//
// --------------------------------------------------------------------
@#for name in [ "c", "w", "i", "p", "R", "z", "s", "a", "me", "m", "t" ]
    fi@{name} = rho@{name}*fi@{name}(-1) - (sig@{name}/100)*eps@{name};
@#endfor
// --------------------------------------------------------------------
//
// D - Observables (linked only to IM equations)
//
// --------------------------------------------------------------------

end;
// --------------------------------------------------------------------
// E N D  O F  M O D E L  B L O C K
// --------------------------------------------------------------------


// --------------------------------------------------------------------
// S T E A D Y  S T A T E  B L O C K
// --------------------------------------------------------------------
//
// Step #1 : define local values
//
// --------------------------------------------------------------------
steady_state_model;
// --------------------------------------------------------------------
// IM Model -> + steady-state restrictions
// --------------------------------------------------------------------
// --------------------------------------------------------------------
p0_IM         = 1;
muz0          = 0;
fss0_IM       = cstf/100;
sss0_IM       = csts/100;
Rss0          = cstR/100;
piss0         = 0;
pi_star0      = 0;
Delta0_IM     = 1;
s0_IM         = csts/100;
f0_IM         = cstf/100;
r0_IM         = (1+Rss0)/(1+piss0)-1;
betaF0        = exp(sig*muz0)*(1+piss0)/(1+Rss0);
mu0           = theta/(theta-1);
MF0_IM        = betaF0/exp(sig*muz0);
rho0          = s0_IM/(1-f0_IM);
n0_IM         = f0_IM/(s0_IM+f0_IM);
m0_IM         = rho0*n0_IM;
neu0_IM       = s0_IM*n0_IM;
nuu0_IM       = 1-n0_IM-neu0_IM;
tau0          = (1-n0_IM)/(n0_IM/reprat); 
v0_IM         = m0_IM^(1/(1-chi))/(1-(1-rho0)*n0_IM)^(chi/(1-chi));
lambda0_IM    = m0_IM/v0_IM;
bF_bW0        = s0_IM*(ceu_ce)^(-sig) + (1-s0_IM);
betaW0        = betaF0/bF_bW0;
psi0          = 1; 
alpo0         = (Omega + (psi0)*(1-Omega));
y0_IM         = exp(-muz0)*(alpo0)*n0_IM;
kappav0       = kv_y*y0_IM/v0_IM;
mc0_IM        = 1/mu0;
q0_IM         = mc0_IM;
w0_IM         = q0_IM -kappav0/lambda0_IM*(1-(1-rho0)*MF0_IM); 
wbar0         = w0_IM*(exp(-muz0)/(1+piss0))^(gammaw/(gammaw-1));
bu0_IM        = reprat*w0_IM;
cuu0_IM       = bu0_IM;
ceu0_IM       = ((1-tau0)*w0_IM-(1-s0_IM-exp(muz0)/(1+r0_IM))*bu0_IM)/(1/ceu_ce - (1-s0_IM-exp(muz0)/(1+r0_IM)) );
lameu0_IM     = ( (1/exp(muz0))*ceu0_IM )^(-sig) ;
ce0_IM        = ceu0_IM/ceu_ce;
c0_IM         = y0_IM - kappav0*v0_IM;
cF0_IM        = (c0_IM - Omega*( n0_IM*ce0_IM + (1-f0_IM)*(1-n0_IM)*cuu0_IM + s0_IM*n0_IM*ceu0_IM ))/(1-Omega);
lamF0_IM      = ( (1/exp(muz0))*cF0_IM )^(-sig) ;
x10_IM        = lamF0_IM*mc0_IM*y0_IM/(1-alpha*betaF0*(1+piss0)^theta);
x20_IM        = lamF0_IM*y0_IM/(1-alpha*betaF0*(1+piss0)^(theta-1));
lame0_IM      = ( (1/exp(muz0))*ce0_IM )^(-sig);
MW0_IM        = betaW0*exp(-sig*(muz0))*( (1-s0_IM)*lame0_IM + s0_IM*lameu0_IM )/lame0_IM;
ae0_IM        = exp(muz0)/(1+r0_IM)*(ceu0_IM-bu0_IM);
Ae0_IM        = exp(-muz0)*(1-s0_IM)*n0_IM*ae0_IM;
aF0_IM        = -Omega/(1-Omega)*Ae0_IM;

// --------------------------------------------------------------------
// CM Model -> Imposing the parameters obtained above
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// CM Model -> Imposing the parameters obtained above
// --------------------------------------------------------------------
p0_HM         = 1;
ae0_HM        = 0;
Ae0_HM        = 0;
r0_HM         = (1+Rss0)/(1+piss0)-1;
MF0_HM        = betaF0/exp(sig*muz0);
MW0_HM        = betaW0/exp(sig*muz0);
w0_HM         = wbar0*(exp(-muz0)/(1+piss0))^(gammaw/(1-gammaw));
f0_HM         = cstf/100;
s0_HM         = rho0*(1-f0_HM);
n0_HM         = f0_HM/(s0_HM+f0_HM);
m0_HM         = rho0*n0_HM;
v0_HM         = m0_HM^(1/(1-chi))/(1-(1-rho0)*n0_HM)^(chi/(1-chi));
lambda0_HM    = m0_HM/v0_HM;
q0_HM         = kappav0/lambda0_HM*(1-(1-rho0)*MF0_HM)+w0_HM;
y0_HM         = exp(-muz0)*n0_HM;
mc0_HM        = 1/mu0;
neu0_HM       = s0_HM*n0_HM;
nuu0_HM       = 1-n0_HM-neu0_HM;
bu0_HM        = reprat*w0_HM;
cuu0_HM       = w0_HM*n0_HM;
ceu0_HM       = w0_HM*n0_HM;
ce0_HM        = w0_HM*n0_HM;
lame0_HM      = ( (1/exp(muz0))*ce0_HM )^(-sig) ;
lameu0_HM     = ( (1/exp(muz0))*ceu0_HM )^(-sig) ;
lamuu0_HM     = ( (1/exp(muz0))*cuu0_HM )^(-sig) ;
c0_HM         = y0_HM - kappav0*v0_HM;
cF0_HM        = (c0_HM - Omega*( n0_HM*ce0_HM + (1-f0_HM)*(1-n0_HM)*cuu0_HM + s0_HM*n0_HM*ceu0_HM ))/(1-Omega);
lamF0_HM      = ( (1/exp(muz0))*cF0_HM )^(-sig) ;
x10_HM        = lamF0_HM*mc0_HM*y0_HM/(1-alpha*betaF0*(1+piss0)^theta);
x20_HM        = lamF0_HM*y0_HM/(1-alpha*betaF0*(1+piss0)^(theta-1));
aF0_HM        = 0;
// --------------------------------------------------------------------
//
// Step #2 : assign these local values
//
// --------------------------------------------------------------------
@#for name in [ "IM","HM"]
    r_@{name}         = r0_@{name};
    R_@{name}         = Rss0;
    MF_@{name}        = MF0_@{name};
    MW_@{name}        = MW0_@{name};
    lamF_@{name}      = lamF0_@{name};
    lame_@{name}      = lame0_@{name};
    lameu_@{name}     = lameu0_@{name};
    cF_@{name}        = cF0_@{name};
    ce_@{name}        = ce0_@{name};
    ceu_@{name}       = ceu0_@{name};
    cuu_@{name}       = cuu0_@{name};
    aF_@{name}        = aF0_@{name};
    ae_@{name}        = ae0_@{name};
    Ae_@{name}        = Ae0_@{name};
    lambda_@{name}    = lambda0_@{name};
    m_@{name}         = m0_@{name};
    v_@{name}         = v0_@{name};
    s_@{name}         = s0_@{name};
    f_@{name}         = f0_@{name};
    n_@{name}         = n0_@{name};
    neu_@{name}       = neu0_@{name};
    nuu_@{name}       = nuu0_@{name};
    rho_@{name}       = rho0;
    y_@{name}         = y0_@{name};
    mc_@{name}        = mc0_@{name};
    w_@{name}         = w0_@{name};
    q_@{name}         = q0_@{name};
    bu_@{name}        = bu0_@{name};
    x1_@{name}        = x10_@{name};
    x2_@{name}        = x20_@{name};
    pi_star_@{name}   = pi_star0;
    Delta_@{name}     = Delta0_IM;
    pi_@{name}        = piss0;
    c_@{name}         = c0_@{name};
    p_@{name}         = p0_@{name};
@#endfor

// --------------------------------------------------------------------
// Shocks
// --------------------------------------------------------------------
@#for name in [ "c", "w", "i", "p", "R", "z", "s", "a", "me", "m", "t" ]
    fi@{name} = 0;
@#endfor
end;
// --------------------------------------------------------------------
// E N D  O F  S T E A D Y  S T A T E  B L O C K
// --------------------------------------------------------------------

// --------------------------------------------------------------------
// WARNING: if you want to check the steady state you can uncomment
// the lines below. Please, re-comment them once you're convinced
// that everything works fine
// --------------------------------------------------------------------
//options_.dynatol.f = 1e-10;
//steady;
//check;
//return;

// --------------------------------------------------------------------
// S H O C K  B L O C K
// --------------------------------------------------------------------
shocks;
var  epsp; stderr 1;
end;
// --------------------------------------------------------------------
// E N D  O F  S H O C K  B L O C K
// --------------------------------------------------------------------
