// --------------------------------------------------------------------
//  V A R I A B L E  D E C L A R A T I O N
// --------------------------------------------------------------------
var
// --------------------------------------------------------------------
// model variables
// --------------------------------------------------------------------
@#for name in ["IM","HM"]
    r_@{name} R_@{name} MF_@{name} MW_@{name} lamF_@{name} lame_@{name} 
    lameu_@{name} cF_@{name} ce_@{name} ceu_@{name} cuu_@{name} c_@{name}
    aF_@{name} ae_@{name} Ae_@{name} lambda_@{name} m_@{name} v_@{name} 
    f_@{name} s_@{name} n_@{name} neu_@{name} nuu_@{name} rho_@{name}
    y_@{name} w_@{name} q_@{name} bu_@{name} Delta_@{name} pi_@{name} 
    pi_star_@{name} x1_@{name} x2_@{name} mc_@{name} p_@{name}
@#endfor
// --------------------------------------------------------------------
// shocks
// --------------------------------------------------------------------
@#for name in [ "c", "w", "i", "p", "R", "z", "s", "a", "me", "m", "t" ]
    fi@{name}
@#endfor
// --------------------------------------------------------------------
// observables
// --------------------------------------------------------------------

;

// --------------------------------------------------------------------
//  E X O G E N O U S  V A R I A B L E  D E C L A R A T I O N
// --------------------------------------------------------------------
varexo
@#for name in [ "c", "w", "i", "p", "R", "z", "s", "a", "me", "m", "t" ]
    eps@{name} $\epsilon_{@{name}}$
@#endfor
;

// --------------------------------------------------------------------
//  P A R A M E T E R S  D E C L A R A T I O N
// --------------------------------------------------------------------
parameters
// --------------------------------------------------------------------
// structural parameters
// --------------------------------------------------------------------
sig $\sigma$ delta $\delta$ theta $\theta$ chi $\chi$ hF $h$ Omega $\Omega$
kv_y ceu_ce
reprat upsilon UPSILON phi
nu $\nu_u$ ni $\nu_i$ alpha $\alpha$ gammap $\gamma_p$ 
gammaw $\gamma_w$ psin $\psi_n$
rhoTR $\rho$ api $a_\pi$ ay $a_y$
cstP cstW cstY cstR cstf csts cstdf cstds csts60 lambdaHP
// --------------------------------------------------------------------
// shocks
// --------------------------------------------------------------------
@#for name in [ "c", "w", "i", "p", "R", "z", "s", "a", "me", "m", "t" ]
    rho@{name} $\rho_{@{name}}$ sig@{name} $\sigma_{@{name}}$
@#endfor
;

