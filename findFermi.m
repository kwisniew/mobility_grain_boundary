function Ef=findFermi(Na)
    
    q=1.6021766208*10^(-19);
    T=300;%K
    eps = 13*8.85*10^-12;%A^2 s^4 / kh m^3 
    epsR = 12;
    L=1e-6;%m
    Vb = q*Na*L^2/(8*eps*epsR);
    Ef = fzero(@charge_neut_cond,-q*Vb);
    
    function cnc_expression = charge_neut_cond(Ef)
      
        k = 1.38064852*10^(-23); %J/K
        e_eff_mass = 0.09;
        h_eff_mass = 0.72;
        Nc = 2.5*10^19*((e_eff_mass)^(3/2))*(T/300)^(3/2);
        Nv= 2.5*10^19*((h_eff_mass)^(3/2))*(T/300)^(3/2);
        Eg = 1*q;
        ni = sqrt(Nv*Nc)*exp(-Eg/(2*k*T))*1e6;

        %szerokoœæ ziarna:
        L=1e-6;%m
        %szerokoœc granicy miêdzy ziarnami:
        delta=1e-7;%m  
        deltaE = 0.05*q;

        W = sqrt(2*eps*epsR*Vb/(q*Na));
        Et=0.01*q;
        delta_Et = 0.08*q;
        delta_eps = delta_Et/(2*sqrt( log(2) ));
        Qt=1e12*1e4;

        
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
        
%         if W>L
        cnc_expression = Na*L -  Qt*(k*T/(sqrt(pi)*delta_eps))*(Q1-Q2) - delta*ni*exp(-q*Ef/(k*T))*exp(-(q*Vb-deltaE)/(k*T));
%         else
%             cnc_expression = 2*W*Na -  Qt*(k*T/(sqrt(pi)*delta_eps))*(Q1-Q2) - Na*exp(-(q*Vb-deltaE)/(k*T));
%         end
    end
end