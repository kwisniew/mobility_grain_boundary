%PARAMETRY
    %Sta³e Fizyczne
    %sta³a boltzmana
    k = 1.38064852*10^(-23); %J/K
    %przenikalnoœæ elektryczna
    eps = 13*8.85*10^-12;%A^2 s^4 / kh m^3 
    % ³adunek elementarny
    q=1.6021766208*10^(-19); %C

    
    %parametry mikroskopowe
    % sta³a sieciowa, potrzebne do "extendent state mobility"
    a = 2e-10; %m
    % jump frequency, potrzebne do "extendent state mobility"
    ni_jump = 1e15; %s^-1
    % phonon frrequency, potrzebne do "hopping mobility"
    ni_ph = 1e13; %s^-1
    % hopping distance, potrzebne do "hopping mobility" 
    R = 2e-10; %m
    % deltaE - band tail? Davis-Mott energy band? Table II, potrzebne do "grain boundary mobility ugb"
    deltaE = 0.05*q; %J
    % mobility shoulder? Table II, potrzebne do "grain boundary mobility ugb"
    deltaE_prim = 0.07*q; %J
    
    
    %parametry pó³przewodnika
    %ruchliwoœæ w krysztale
    uc = 0.02*1e-4; %m^2/V*s
    % domieszkoanie
    Na=1e21; %m^-3
    % gêstoœæ dziur w wartwie neutralnej p(L/2)
    p_L2 =Na;
    % przenikalnoœæ elektryczna wzglêdna
    epsR = 12;
    %masa efektywna elektronów
    e_eff_mass = 0.09;  
    %masa efektywna dziur
    h_eff_mass = 0.72;
    %efektywna gêstoœæ stanów pas przewodnictwa
    Nc = 2.5*10^19*((e_eff_mass)^(3/2))*(T/300)^(3/2);
    %efektywna gêstoœæ stanów pas walencyjny
    Nv= 2.5*10^19*((h_eff_mass)^(3/2))*(T/300)^(3/2);
    %przerwa energetyczna
    Eg = 1*q;
    % gêstoœæ noœników w przewodnik samoistnym
    ni = sqrt(Nv*Nc)*exp(-Eg/(2*k*T))*1e6;

    
    %parametry ogólne
    %temperatura
    T=300;%K
    %szerokoœæ ziarna
    L=3e-6;%m
    %szerokoœæ granicy ziarna
    delta=1e-7;%m
    
    
    %parametry pomocnicze do granicy ziaren
    % kappa, Table I, potrzebne do gammy
    kappa = sqrt( (deltaE+deltaE_prim)/(k*T) );
    % gamma (Appendix)
    gamma = erfc(kappa)+2*kappa*exp(-kappa^2)/sqrt(pi);
    % preexponential mobility, równanie (7)
    u0 = ni_ph*q*R^2/(6*k*T); %m/Vs

        
    %parametry granicy ziaren
    %bariera na granicy ziarna
    Vb=0.1; %J = 0.1eV
    % gêstoœæ dziur na granicy ziarna, Appendix
    pgb = p_L2*exp(-(q*Vb - deltaE)/(k*T));%m^-3
    %extendent state mobility
    uext = q*a^2*ni_jump/(6*k*T); % m/Vs
    % hopping barrier H(E)
    H = k*T; %J
    %hopping mobility    
    uhop = u0*exp(-(H/(k*T))); %m/Vs
    % ruchliwoœæ dziór na granicy ziarna, Appendix
    ugb = gamma*uext+(1-gamma)*uhop;

  
    %szerokoœæ warstwy zubo¿onej granicy ziaren
    W = sqrt(2*eps*epsR*Vb/(q*Na));
    
    %Fgb
    Fgb = (delta/L)*(uc/ugb)*(p_L2/pgb);
    
    %Fc0
    Fc = 1- 2*W/L + ( (2*W-delta)/L )* exp( q*Vb/(k*T) )*dawson( sqrt(q*Vb/(k*T)) )/ sqrt(q*Vb/(k*T));
    
    %ruchliwoœæ z lotu ptaka
    u = uc/(Fgb + Fc);
    
    %chcemy zobaczyæ zmiennoœc wzglêdem "Na", "Vb" oraz "delta"
    % Od Vb zale¿y "pgb", "W" oraz "Fc"
    % Od Na zale¿y jedynie "W" --> uznajê, ¿e "p" równie¿ równa siê "Na"
    % Od delty zale¿y "Fgb" oraz "Fc"
  
    fprintf('Gêstoœæ noœników na granicy ziarna:  %e cm^-3\n', pgb_fun(p_L2,Vb,deltaE,T)/1e6);
    fprintf('Warstwa zubo¿ona:  %f um\n', W_fun(eps,epsR,Vb,Na)*1e6);
    fprintf('Fc: %f\n', Fc_fun(W,L,delta,Vb,T));
    fprintf('Fgb: %f\n', Fgb_fun(delta,L,uc,ugb,p_L2,pgb));
    fprintf('u = %f cm^2/Vs\n', u_fun(uc,Fgb,Fc)*1e4);

    
%ZNALEZIENIE Ef NUMERYCZNIE    
    %zmienne do obliczenia przy danym Na:
    % p_L2, pgb, Ef (fermi level), Vb, delta_Qt oraz W (dla przypadku, gdy ziarno ju¿ nie jest ca³kowicie zubo¿one)

    % obliczam dla jakiej maksymalnej koncentracji akceptorów pó³przewodnik nie jest jeszcze zdegenerowany
    % czyli gdy Ef = Eg/2 (póŸniej obliczenia nie maj¹ sensu)
    max_non_degenerate_Na = ni*exp(Eg/(2*k*T));

    Na   = exp(20*log(10):0.1:log(max_non_degenerate_Na));
    pgb  = zeros(length(Na),1)';
    p_L2 = zeros(length(Na),1)';
    W    = zeros(length(Na),1)';
    Vb   = zeros(length(Na),1)';
    Ef   = zeros(length(Na),1)';
    Fermi_graniczny = k*T*log(ni./Na)/q;
    
    %sklejanie gdy idziemy od ma³ego domieszkowania w górê
    %sklejamy poziomem fermiego
%     i=1;
%     is_fermi_small_enough = true;
%     while i <= length(Na) && is_fermi_small_enough == true
%         Ef(i)   = findFermi(Na(i));
%         p_L2(i) = p_L2_fun(ni,Ef(i),T);
%         W(i)    = L/2;
%         Vb(i)   = Vb_fun(Na(i),W(i));
%         pgb(i)  = pgb_fun(p_L2(i),Vb(i),deltaE,T);
%       
%         if Ef(i) <= Fermi_graniczny(i) || -Ef(i) > 0.5*Eg/q 
%             is_fermi_small_enough = false;
%         end
%         i=i+1;
%             
%     end
%     if is_fermi_small_enough == false
%         i=i-1;
%         is_fermi_small_enough = true;
%         while i <= length(Na) && is_fermi_small_enough == true
%             Ef(i)   = Fermi_graniczny(i);
%             p_L2(i) = p_L2_fun(ni,Ef(i),T);
%             W(i)    = findW(Na(i));
%             Vb(i)   = Vb_fun(Na(i),W(i));
%             pgb(i)  = pgb_fun(p_L2(i),Vb(i),deltaE,T);
% 
%             if -Ef(i) > 0.5*Eg/q 
%                 is_fermi_small_enough = false;
%             end
%             i=i+1;
% 
%         end
%     end

%sklejanie gdy idziemy od du¿ych gêstoœci domieszkowania
%sklejanie wielkoœci¹ warstwy zubo¿onej
i=length(Na);
is_W_small_enough = true;
while i > 0 && is_W_small_enough == true
    
    Ef(i)   = Fermi_graniczny(i);
    p_L2(i) = p_L2_fun(ni,Ef(i),T);
    W(i)    = findW(Na(i));
    Vb(i)   = Vb_fun(Na(i),W(i));
    pgb(i)  = pgb_fun(p_L2(i),Vb(i),deltaE,T);

    if W(i) > L/2 
        is_W_small_enough = false;
    end
    i=i-1;
    
end

if is_W_small_enough == false
    i=i+1;
    while i > 0
        Ef(i)   = findFermi(Na(i));
        p_L2(i) = p_L2_fun(ni,Ef(i),T);
        W(i)    = L/2;
        Vb(i)   = Vb_fun(Na(i),W(i));
        pgb(i)  = pgb_fun(p_L2(i),Vb(i),deltaE,T);

        i=i-1;

    end
end
    
    
    
%WYKRESY:
subplot(1,3,1)
loglog(Na/1e6,u_fun(uc, ...
                Fgb_fun(delta,L,uc,ugb,p_L2,pgb), ...
                Fc_fun(W,L,delta,Vb,T))*1e4)
title('ruchliwoœæ ~ domieszkowania')

subplot(1,3,2)
loglog(Na/1e6,W)
title('Warstwa zubo¿ona ~ domieszkowania')

subplot(1,3,3)
loglog(Na/1e6,Ef)
title('Poziom Fermiego ~ domieszkowania')
                
% u_fun(uc, ...
%       Fgb_fun(delta,L,uc,ugb,p_L2(1),pgb(1)), ...
%       Fc_fun(W(1),L,delta,Vb(1),T))*1e4  


%     Vb_x = [0.01:0.01:0.35];
%     tmp = u_fun(uc,Fgb,Fc_fun(W_fun(eps,epsR,Vb_x,Na),L,delta,Vb,T))*1e4;
%     %plot(Vb_x,   u_fun(uc,Fgb,Fc_fun(W_fun(eps,epsR,Vb_x,Na),L,delta,Vb,T))*1e4   )
%     plot(Vb_x,u_fun(uc, ...
%                     Fgb_fun(delta,L,uc,ugb,p_L2,pgb_fun(p_L2,Vb_x,deltaE,T)), ...
%                     Fc_fun(W_fun(eps,epsR,Vb_x,Na),L,delta,Vb_x,T))*1e4)
                
   %figure
   %semilogy(Vb_x,u_fun(uc, ...
   %                    Fgb_fun(delta,L,uc,ugb,p,pgb_fun(p,Vb_x,deltaE,T)), ...
   %                    Fc_fun(W_fun(eps,epsR,Vb_x,Na),L,delta,Vb_x,T))*1e4)
            
%     Na_x = [1e21:1e21:1e23];
%     plot(Na_x,u_fun(uc, ...
%                     Fgb_fun(delta,L,uc,ugb,Na_x,pgb_fun(Na_x,Vb,deltaE,T)), ...
%                     Fc_fun(W_fun(eps,epsR,Vb,Na_x),L,delta,Vb,T))*1e4)
%                 
%     plot(Na_x,pgb_fun(Na_x,Vb,deltaE,T))
%     plot(Na_x,Fc_fun(W_fun(eps,epsR,Vb,Na_x),L,delta,Vb,T))
%     plot(Na_x,Fgb_fun(delta,L,uc,ugb,Na_x,pgb_fun(Na_x,Vb,deltaE,T)))           

    
%FUNKCJE
    function Vb = Vb_fun(Na,W)
        q    = 1.6021766208*10^(-19);
        eps  = 13*8.85*10^-12;%A^2 s^4 / kh m^3 
        epsR = 12;
        Vb   = q*Na*W^2/(2*eps*epsR);
    end  
    function p_L2 = p_L2_fun(ni,Ef,T)
        q    = 1.6021766208*10^(-19);
        k = 1.38064852*10^(-23); %J/K
        p_L2 = ni*exp(-q*Ef/(k*T));
    end    
    function pgb = pgb_fun(p_L2,Vb,deltaE,T)
        q = 1.6021766208*10^(-19);
        k = 1.38064852*10^(-23); %J/K
        pgb = p_L2.*exp(-(q*Vb - deltaE)/(k*T));
    end    
    
    function W = W_fun(eps,epsR,Vb,Na)
        q=1.6021766208*10^(-19);
        W = sqrt(2*eps*epsR*Vb.*(q*Na).^-1);
    end
    
   function Fc = Fc_fun(W,L,delta,Vb,T)
        q=1.6021766208*10^(-19);
        k = 1.38064852*10^(-23); %J/K
        Fc = 1- 2*W/L + ( (2*W-delta)/L ).* exp( q*Vb/(k*T) ).*dawson( sqrt(q*Vb/(k*T)) ).* sqrt(q*Vb/(k*T)).^-1;
   end
    
   function Fgb = Fgb_fun(delta,L,uc,ugb,p_L2,pgb)
        Fgb = (delta/L)*(uc/ugb)*(p_L2.*(pgb.^-1));
   end
   
   function u = u_fun(uc,Fgb,Fc)
        u = uc *(Fgb + Fc).^-1;
   end

    
