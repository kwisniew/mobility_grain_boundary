function W=findW(Na, T, delta, L, Qt, Et, delta_Et, Eg, epsR)

    W = fzero(@charge_neut_cond,L/2);
    
    function cnc_expression = charge_neut_cond(W)
        
        q=1.6021766208*10^(-19);
        eps = 8.85*10^-12;%A^2 s^4 / kg m^3 
        k = 1.38064852*10^(-23); %J/K
        e_eff_mass = 0.09;
        h_eff_mass = 0.72;
        Nc = 2.5*10^19*((e_eff_mass)^(3/2))*(T/300)^(3/2);
        Nv= 2.5*10^19*((h_eff_mass)^(3/2))*(T/300)^(3/2);
        ni = sqrt(Nv*Nc)*exp(-Eg/(2*k*T))*1e6;
        Ef = k*T*log(ni./Na)/q;

        deltaE = 0.05*q;
        delta_eps = delta_Et/(2*sqrt( log(2) ));
        Vb = q*Na*W^2/(2*eps*epsR);
        
        
        b1 = (q*Ef+Et+q*Vb)/(k*T);
        Q1 =  Q_integral_fun(b1,5.12634,2.06901,0.11965,0.0230097,0.721425,0.0829264,0.172877);
         
        b2 = (Et)/(k*T);
        Q2 =  Q_integral_fun(b2,5.12634,2.06901,0.11965,0.0230097,0.721425,0.0829264,0.172877);
       
        cnc_expression = 2*W*Na -  Qt*(k*T/(sqrt(pi)*delta_eps))*(Q1-Q2) - delta*ni*exp(-q*Ef/(k*T))*exp(-(q*Vb-deltaE)/(k*T));
        
    end
end

function Q_integral = Q_integral_fun(b,A,B,C,D,E,F,G)
    Q_integral =  A-B*exp((C-D*b)*b)+exp((-C-D*b)*b)*(B-B*erf(E-F*b)) + ...
                    B*exp((C-D*b)*b)*erf(E+F*b)-A*erf(G*b);
end
%  ca³ka dla T=300, delta_Et=0.08q
%  Q_integral_fun(b2,1.64696,0.844494,0.386513,0.0743295,0.416573,0.46392,0.5381);
% 
%         Q1 =  1.64696 - 0.844494 * exp((0.386513 - 0.0743295 * b1)* b1) + ... 
%               exp((-0.386513 - 0.0743295 *b1) *b1) * (0.844494 - 0.844494 * erf(0.416573 - 0.46392*b1)) + ...
%               0.844494 * exp((0.386513 - 0.0743295 *b1) *b1)* erf(0.416573 + 0.46392 * b1) - ... 
%               1.64696 * erf(0.5381 * b1);
%         Q2 =  1.64696 - 0.844494 * exp((0.386513 - 0.0743295 * b2)* b2) + ... 
%               exp((-0.386513 - 0.0743295 *b2) *b2) * (0.844494 - 0.844494 * erf(0.416573 - 0.46392*b2)) + ...
%               0.844494 * exp((0.386513 - 0.0743295 *b2) *b2)* erf(0.416573 + 0.46392 * b2) - ... 
%               1.64696 * erf(0.5381 * b2);



%a = k*T/delta_eps;
%Q1 =  (0.443113 * exp((0.0676+(a^2)*(-0.52-0.1*b1)*b1)/(0.1+a^2))*(1+erf((-0.26+(a^2)*b1)/sqrt(0.1+a^2))))/sqrt(0.1+a^2)+exp(-a^2*b1^2)*(-((exp(a^2*b1^2)* sqrt(pi)*(-1+erf(sqrt(a^2)*b1)))/(2*sqrt(a^2)))+(0.443113*exp((0.26+(a^2)*b1)^2/(0.1+a^2))*(-1.0+erf((0.26+(a^2)*b1)/sqrt(0.1+a^2))))/sqrt(0.1+a^2));
%         Q2 =  (0.443113 * exp((0.0676 + (a^2)*(-0.52-0.1*b2)*b2)/(0.1+a^2))* ...
%               (1 + erf((-0.26 + (a^2)*b2)/sqrt(0.1 + a^2)))) / sqrt(0.1 + a^2) + ...
%                exp(-a^2*b2^2)* ...
%                 (-((exp(a^2*b2^2)* sqrt(pi)*(-1 + erf(sqrt(a^2)*b2)))/(2*sqrt(a^2))) + ...
%                 (0.443113*exp((0.26 + (a^2)*b2)^2/(0.1 + a^2))*(-1. + erf((0.26 + (a^2)*b2)/sqrt(0.1 + a^2))))/...
%                    sqrt(0.1 + a^2));        