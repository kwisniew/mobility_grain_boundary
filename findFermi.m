function Ef=findFermi(Na, T, delta, L, Qt, Et, delta_Et, Eg, epsR)
    
    q=1.6021766208*10^(-19);
    eps = 8.85*10^-12;%A^2 s^4 / kg m^3 
    Vb = q*Na*L^2/(8*eps*epsR);
    %options = optimset('Display','iter');
    Ef = fzero(@charge_neut_cond,-Vb);%,options);
    
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
        switch T
            case 100
                Q1 = Q_integral_fun(b1,5.12634,2.06901,0.11965,0.0230097,0.721425,0.0829264,0.172877);
            case 125
                Q1 = Q_integral_fun(b1,4.10107,1.83414,0.165529,0.0318326,0.721425,0.0829264,0.172877);;
            case 150
                Q1 = Q_integral_fun(b1,3.41756,1.62323,0.209078,0.0402073,0.635766,0.16443,0.259316);
            case 175
                Q1 = Q_integral_fun(b1,2.92934,1.44107,0.248498,0.0477881,0.594098,0.209139,0.302535);
            case 200
                Q1 = Q_integral_fun(b1,2.56317,1.28669,0.283147,0.0544514,0.554895,0.255136,0.345754);
            case 250
                Q1 = Q_integral_fun(b1,2.05054,1.04737,0.338683,0.0651313,0.485502,0.348797,0.432193);
            case 300
                Q1 = Q_integral_fun(b1,1.64696,0.844494,0.386513,0.0743295,0.416573,0.46392,0.5381);
            otherwise
                fprintf('CALKA DLA TEJ TEMPERATURY NIE ZAIMPELEMNTOWANA\n');
        end
        
        

        b2 = (Et)/(k*T);
        switch T
            case 100
                Q2 = Q_integral_fun(b2,5.12634,2.06901,0.11965,0.0230097,0.721425,0.0829264,0.172877);
            case 125
                Q2 = Q_integral_fun(b2,4.10107,1.83414,0.165529,0.0318326,0.721425,0.0829264,0.172877);;
            case 150
                Q2 = Q_integral_fun(b2,3.41756,1.62323,0.209078,0.0402073,0.678831,0.121922,0.216096);
            case 175
                Q2 = Q_integral_fun(b2,2.92934,1.44107,0.248498,0.0477881,0.594098,0.209139,0.302535);
            case 200
                Q2 = Q_integral_fun(b2,2.56317,1.28669,0.283147,0.0544514,0.554895,0.255136,0.345754);
            case 250
                Q2 = Q_integral_fun(b2,2.05054,1.04737,0.338683,0.0651313,0.485502,0.348797,0.432193);
            case 300
                Q2 = Q_integral_fun(b2,1.64696,0.844494,0.386513,0.0743295,0.416573,0.46392,0.5381);
        end
        Q_dep_reg = Na*L;
        Q_trap    = Qt*(k*T/(sqrt(pi)*delta_eps))*(Q1-Q2);
        Q_p_gb    = delta*ni*exp(-q*Ef/(k*T))*exp(-(q*Vb-deltaE)/(k*T));
        %max_Q     = max([Q_dep_reg,Q_trap,Q_p_gb]);
        cnc_expression = (Q_dep_reg - Q_trap - Q_p_gb);
        
    end
end

%heurystyka przepisywania parametrów A,B,...:
%najpierw przepisz 4 pierwsze parametry jak leci, póŸniej dwa parametry z
%pierwszego erfa, a póŸniej parametr z ostatniego erfa
function Q_integral = Q_integral_fun(b,A,B,C,D,E,F,G)
    Q_integral =  A-B*exp((C-D*b)*b)+exp((-C-D*b)*b)*(B-B*erf(E-F*b)) + ...
                    B*exp((C-D*b)*b)*erf(E+F*b)-A*erf(G*b);
end

%ca³ka dla T=100, delta_Et=0.08q
%Q_integral_fun(b1,5.12634,2.06901,0.11965,0.0230097,0.721425,0.0829264,0.172877);

%ca³ka dla T=150, delta_Et=0.08q
%Q_integral_fun(b2,3.41756,1.62323,0.209078,0.0402073,0.635766,0.16443,0.259316);

%ca³ka dla T=200, delta_Et=0.08q
%Q_integral_fun(b2,2.56317,1.28669,0.283147,0.0544514,0.554895,0.255136,0.345754)

% ca³ka dla T=250, delta_Et=0.08q
% Q_integral_fun(b1,2.05054,1.04737,0.338683,0.0651313,0.485502,0.348797,0.432193);


%  ca³ka dla T=300, delta_Et=0.08q
%
%         Q1 =  1.64696 - 0.844494 * exp((0.386513 - 0.0743295 * b1)* b1) + ... 
%               exp((-0.386513 - 0.0743295 *b1) *b1) * (0.844494 - 0.844494 * erf(0.416573 - 0.46392*b1)) + ...
%               0.844494 * exp((0.386513 - 0.0743295 *b1) *b1)* erf(0.416573 + 0.46392 * b1) - ... 
%               1.64696 * erf(0.5381 * b1);
%         Q2 =  1.64696 - 0.844494 * exp((0.386513 - 0.0743295 * b2)* b2) + ... 
%               exp((-0.386513 - 0.0743295 *b2) *b2) * (0.844494 - 0.844494 * erf(0.416573 - 0.46392*b2)) + ...
%               0.844494 * exp((0.386513 - 0.0743295 *b2) *b2)* erf(0.416573 + 0.46392 * b2) - ... 
%               1.64696 * erf(0.5381 * b2);


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