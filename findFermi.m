function Ef=findFermi(Na, T, delta, L, Qt, Et, delta_Et, Eg, epsR)
    
    q=1.6021766208*10^(-19);
    eps = 8.85*10^-12;%A^2 s^4 / kg m^3 
    Vb = q*Na*L^2/(8*eps*epsR);
    %options = optimset('Display','iter');
    Ef = fzero(@charge_neut_cond,0.0);%,options);
    
    function cnc_expression = charge_neut_cond(Ef)
      
        k = 1.38064852*10^(-23); %J/K
        e_eff_mass = 0.09;
        h_eff_mass = 0.72;
        Nc = 2.5*10^19*((e_eff_mass)^(3/2))*(T/300)^(3/2);
        Nv= 2.5*10^19*((h_eff_mass)^(3/2))*(T/300)^(3/2);
        ni = sqrt(Nv*Nc)*exp(-Eg/(2*k*T))*1e6;

        deltaE = 0.05*q;
        delta_eps = delta_Et/(2*sqrt( log(2) ));

        
        b1 = (q*Ef+Et+q*Vb)/(k*T);
        
        Q1 =  1.64696 - 0.844494 * exp((0.386513 - 0.0743295 * b1)* b1) + ... 
              exp((-0.386513 - 0.0743295 *b1) *b1) * (0.844494 - 0.844494 * erf(0.416573 - 0.46392*b1)) + ...
              0.844494 * exp((0.386513 - 0.0743295 *b1) *b1)* erf(0.416573 + 0.46392 * b1) - ... 
              1.64696 * erf(0.5381 * b1);

        b2 = (Et)/(k*T);                                       
        Q2 =  1.64696 - 0.844494 * exp((0.386513 - 0.0743295 * b2)* b2) + ... 
              exp((-0.386513 - 0.0743295 *b2) *b2) * (0.844494 - 0.844494 * erf(0.416573 - 0.46392*b2)) + ...
              0.844494 * exp((0.386513 - 0.0743295 *b2) *b2)* erf(0.416573 + 0.46392 * b2) - ... 
              1.64696 * erf(0.5381 * b2);

        Q_dep_reg = Na*L;
        Q_trap    = Qt*(k*T/(sqrt(pi)*delta_eps))*(Q1-Q2);
        Q_p_gb    = delta*ni*exp(-q*Ef/(k*T))*exp(-(q*Vb-deltaE)/(k*T));
        %max_Q     = max([Q_dep_reg,Q_trap,Q_p_gb]);
        cnc_expression = (Q_dep_reg - Q_trap - Q_p_gb);
        
    end
end

        %a=0.5186;%k*T/delta_eps;
%         Q1 =  (0.443113 * exp((0.0676 + (a^2)*(-0.52-0.1*b1)*b1)/(0.1+a^2))* ...
%               (1 + erf((-0.26 + (a^2)*b1)/sqrt(0.1 + a^2)))) / sqrt(0.1 + a^2) + ...
%                exp(-a^2*b1^2)* ...
%                 (-((exp(a^2*b1^2)* sqrt(pi)*(-1 + erf(sqrt(a^2)*b1)))/(2*sqrt(a^2))) + ...
%                 (0.443113*exp((0.26 + (a^2)*b1)^2/(0.1 + a^2))*(-1. + erf((0.26 + (a^2)*b1)/sqrt(0.1 + a^2))))/...
%                    sqrt(0.1 + a^2));

%         Q2 =  (0.443113 * exp((0.0676 + (a^2)*(-0.52-0.1*b2)*b2)/(0.1+a^2))* ...
%               (1 + erf((-0.26 + (a^2)*b2)/sqrt(0.1 + a^2)))) / sqrt(0.1 + a^2) + ...
%                exp(-a^2*b2^2)* ...
%                 (-((exp(a^2*b2^2)* sqrt(pi)*(-1 + erf(sqrt(a^2)*b2)))/(2*sqrt(a^2))) + ...
%                 (0.443113*exp((0.26 + (a^2)*b2)^2/(0.1 + a^2))*(-1. + erf((0.26 + (a^2)*b2)/sqrt(0.1 + a^2))))/...
%                    sqrt(0.1 + a^2));