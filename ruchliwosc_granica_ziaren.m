%% PARAMETRY
%TYLKO KILKA PARAMETRÓW WYSTÊPUJE JAKO ZMIENNA W OBLICZENIACH
%S¥ TO: "Na", "Qt", "T", "delta", "Et" i "L"

    % ³adunek elementarny
    q=1.6021766208*10^(-19); %C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %ZMIENNE
    %temperatura (musi byæ od najwiêkszej do najmniejszej)
    % UWAGA: ka¿da zmiana T powoduje, ¿e
    %        trzeba jeszcze raz przeliczyæ ca³ki w Mathematice!!!
%     Temperature = [300,250,200,150,100];
    Temperature = [200,175,150,125,100];
    T=Temperature(5);%K
    %szerokoœæ ziarna
    L=1e-6;%m
    %szerokoœæ granicy ziarna
    delta=2e-9;%m
    % gêstoœæ przestrzenna ³adunku na granicy ziaren 
    Qt=1e10*1e4;%m^-2
    % po³o¿enie Et wzglêdem Ei
    Et=0.3*q;%eV
    % szerokoœæ po³ówkowa 
    % UWAGA: ka¿da zmiana delta_Et (bazowo: 0.083q) powoduje, ¿e
    %        trzeba jeszcze raz przeliczyæ ca³ki w Mathematice!!!
    delta_Et = 0.083*q;
    %przerwa energetyczna
    Eg = 1.20*q;
    % przenikalnoœæ elektryczna wzglêdna
    epsR = 13;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %obliczeniowe
    calculate_thermal_scan = true;
    S=0.025*1e-4;

    %Sta³e Fizyczne    
    %sta³a boltzmana
    k = 1.38064852*10^(-23); %J/K
    %przenikalnoœæ elektryczna
    eps = 8.85*10^-12;%A^2 s^4 / kg m^3 

    
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
    %masa efektywna elektronów
    e_eff_mass = 0.09;  
    %masa efektywna dziur
    h_eff_mass = 0.72;
   
    
    %parametry pomocnicze do granicy ziaren
    % kappa, Table I, potrzebne do gammy
    kappa = sqrt( (deltaE+deltaE_prim)/(k*T) );
    % gamma (Appendix)
    gamma = erfc(kappa)+2*kappa*exp(-kappa^2)/sqrt(pi);
    % preexponential mobility, równanie (7)
    u0 = ni_ph*q*R^2/(6*k*T); %m/Vs

        
    %parametry granicy ziaren
    %extendent state mobility
    uext = 12*1e-4;% q*a^2*ni_jump/(6*k*T); % m^2/Vs 
    % hopping barrier H(E)
    H = k*T; %J
    %hopping mobility    
    uhop = 0.1e-4;%u0*exp(-(H/(k*T))); %m^2/Vs 0.1e-4;%
    % ruchliwoœæ dziur na granicy ziarna, Appendix
    ugb = gamma*uext+(1-gamma)*uhop;
        
    %sprawdzanie
     delta_eps = delta_Et/(2*sqrt( log(2) ));
     deltaE = 0.05*q;  
%% OBLICZENIA NUMERYCZNIE
%zmienne do obliczenia przy danym Na:
% p_L2, pgb, Ef (fermi level), Vb, delta_Qt oraz W (dla przypadku, gdy ziarno ju¿ nie jest ca³kowicie zubo¿one)

% obliczam dla jakiej maksymalnej koncentracji akceptorów pó³przewodnik nie jest jeszcze zdegenerowany
% czyli gdy Ef = Eg/2 (póŸniej obliczenia nie maj¹ sensu) - zob. max_non_degenerate_Na

if calculate_thermal_scan
    number_of_calculation = length(Temperature);
    %UWAGA! poni¿sza wartoœæ bêdzie najmniejsza dla najwiêkszych temperatur
    %       trzeba pamiêtaæ o tym, ¿eby w "Temperature" na pozycji 1 by³a 
    %       najwiêksza temperatura
    ni = ni_fun(Temperature(1),e_eff_mass,h_eff_mass,Eg,k);
    max_non_degenerate_Na = ni*exp(Eg/(2*k*Temperature(1)));

    Na   = exp(19*log(10):0.1:log(max_non_degenerate_Na));
    %Dla du¿ego przybli¿enia na osi Na, mo¿na u¿yæ poni¿szych wartoœci dla
    %akceptorów
%     Na   = 1.3e19*linspace(1,2.5,6000);
%     Na   = 5e19*linspace(1,3,1000);
    pgb  = zeros(length(Na),number_of_calculation)';
    p_L2 = zeros(length(Na),number_of_calculation)';
    W    = zeros(length(Na),number_of_calculation)';
    Vb   = zeros(length(Na),number_of_calculation)';
    Ef   = zeros(length(Na),number_of_calculation)';
    

    %sklejanie gdy idziemy od ma³ego domieszkowania w górê
    %sklejamy poziomem fermiego
    for i=1:number_of_calculation
        ni = ni_fun(Temperature(i),e_eff_mass,h_eff_mass,Eg,k);
        Fermi_graniczny = k*Temperature(i)*log(ni./Na)/q;
        [pgb(i,:),p_L2(i,:),W(i,:),Vb(i,:),Ef(i,:)] = sove_problem_for_all_Na(...
                                    Na,pgb(i,:),p_L2(i,:),W(i,:),Vb(i,:),Ef(i,:),Fermi_graniczny,...
                                    ni,Temperature(i),delta,L,Qt,Et,delta_Et,Eg,epsR,deltaE,q,uc,ugb);
    end
    [Ea, all_mobilities] = calculate_activation_energy(length(Na),Temperature,W,L,delta,Vb,uc,ugb,p_L2,pgb);
else
    number_of_calculation = 1;
    max_non_degenerate_Na = ni*exp(Eg/(2*k*T));

    Na   = exp(19*log(10):0.1:log(max_non_degenerate_Na));
    pgb  = zeros(length(Na),1)';
    p_L2 = zeros(length(Na),1)';
    W    = zeros(length(Na),1)';
    Vb   = zeros(length(Na),1)';
    Ef   = zeros(length(Na),1)';
    Fermi_graniczny = k*T*log(ni./Na)/q;

    %sklejanie gdy idziemy od ma³ego domieszkowania w górê
    %sklejamy poziomem fermiego
    [pgb,p_L2,W,Vb,Ef] = sove_problem_for_all_Na(...
                                    Na,pgb,p_L2,W,Vb,Ef,Fermi_graniczny,...
                                    ni,T,delta,L,Qt,Et,delta_Et,Eg,epsR,deltaE,q,uc,ugb);


end
%% Sklejanie W

%sklejanie gdy idziemy od du¿ych gêstoœci domieszkowania
%sklejanie wielkoœci¹ warstwy zubo¿onej - W
% i=length(Na);
% is_W_small_enough = true;
% while i > 0 && is_W_small_enough == true
%     
%     Ef(i)   = Fermi_graniczny(i);
%     p_L2(i) = p_L2_fun(ni,Ef(i),T);
%     W(i)    = findW(Na(i), T, delta, L, Qt, Et, delta_Et, Eg, epsR);
%     Vb(i)   = Vb_fun(Na(i),W(i));
%     pgb(i)  = pgb_fun(p_L2(i),Vb(i),deltaE,T);
% 
%     if W(i) > L/2 
%         is_W_small_enough = false;
%     end
%     i=i-1;
%     
% end
% 
% if is_W_small_enough == false
%     i=i+1;
%     while i > 0
%         Ef(i)   = findFermi(Na(i), T, delta, L, Qt, Et, delta_Et, Eg, epsR);
%         p_L2(i) = p_L2_fun(ni,Ef(i),T);
%         W(i)    = L/2;
%         Vb(i)   = Vb_fun(Na(i),W(i));
%         pgb(i)  = pgb_fun(p_L2(i),Vb(i),deltaE,T);
% 
%         i=i-1;
% 
%     end
% end
%% WYKRESY:
% subplot(1,3,1)
figure(1)
    semilogx(Na/1e6,all_mobilities*1e4)
    title('ruchliwosc ~ domieszkowania')

figure(2)
    semilogx(Na/1e6,Vb)
    title('bariera ~ domieszkowania')

% loglog(Na/1e6,W)
% title('Warstwa zubo¿ona ~ domieszkowania')
% 
figure(3)
    loglog(Na/1e6,p_L2/1e6)
    title('p(L/2) ~ domieszkowania')

%semilogx(Na/1e6,Ef)
% title('Poziom Fermiego ~ domieszkowania')

if calculate_thermal_scan
%     semilogx(Na(1,:)/1e6,Ea(1,:))
%     xline(8e14);
%     title('Energia aktywacji ~ domieszkowania')
%     my_index=find(abs(Na-8e14*1e6)< 0.2e15*1e6);
%     Na(my_index)/1e6
%     Ea(1,my_index)

%semilogy(1./Temperature, all_mobilities(:,my_index(4))', '*','MarkerSize',10)
%     my_Ea = zeros(2,1);
%     my_Ea(1) = -Ea(1,my_index(4))*q/k;
%     my_Ea(2) =  Ea(2,my_index(4));
%     y_est_mob = polyval(my_Ea,1./Temperature);
% hold on
% semilogy(1./Temperature,exp(y_est),'--','LineWidth',2)

%DANE DOŒWIADCZALNE - po prostu przeklei³em
    amplitudy = [2484, 1282, 618, 316, 164, 77, 39, 20];
    temperatury = [166.08, 148.34, 136.08, 125.03, 120.58, 110.66, 100.73, 90.832];
    c = polyfit(1./temperatury,log(amplitudy),1);
    y_est = polyval(c,1./temperatury);

%     figure
%     yyaxis left
%     semilogy(1./Temperature, all_mobilities(:,my_index(4))', '*',    1./Temperature, exp(y_est_mob),'-')
% 
%     yyaxis right
%     semilogy(1./temperatury,amplitudy, 'x',1./temperatury,exp(y_est),'-')
    % semilogy(1./temperatury,exp(y_est),'--','LineWidth',2)
    % hold off
    
% 3 ro¿ne podejœcia do wyliczania pr¹du nasycenia: 1) zmienne p(L/2) i Vb,
% 2) zmienne tylko Vb (_v2) i 3) zmienna p(L/2), Vb, i ruchliwoœæ
    Ea_Io = zeros(2,length(Na));
    Ea_Io_v2 = zeros(2,length(Na));
    Ea_Io_v3 = zeros(2,length(Na));
    for i=1:length(Na)
        Ea_Io(:,i)    = polyfit(1./Temperature',log( p_L2(:,i).*uc                  .*exp(-q*Vb(:,i)./(k*Temperature')).*Emax(Vb(:,i),0,Na(i),epsR*eps) ),1);
    end
    for i=1:length(Na)
        Ea_Io_v2(:,i) = polyfit(1./Temperature',log( 1e25      *uc*0.0001           .*exp(-q*Vb(:,i)./(k*Temperature')).*Emax(Vb(:,i),0,Na(i),epsR*eps) ),1);
    end
    for i=1:length(Na)
        Ea_Io_v3(:,i) = polyfit(1./Temperature',log(  p_L2(:,i).*all_mobilities(:,i).*exp(-q*Vb(:,i)./(k*Temperature')).*Emax(Vb(:,i),0,Na(i),epsR*eps)),1);
    end

    
    figure(4)
        semilogx(Na/1e6,-Ea_Io_v3(1,:)*k/q)
    %    semilogx(Na/1e6,-Ea_Io_v2(1,:)*k/q)
    %    semilogx(Na/1e6,-Ea_Io(1,:)*k/q)
        yline(0.09);
        title('Ea ~ domieszkowania')
    
    S=0.28*1e-4;%powierzchnia próbki w m2
    %aby mieæ dobre jednoski dla Io/S, w [mA/cm2] muszê pomno¿yæ moje wyniki
    %które równaj¹ siê Io/(qS) przez q, uzyskam jednostki A/m2, muszê
    %jeszcze zamieniæ m2 na cm2 czyli pomno¿yæ wszystko przez 1e-4 i aby
    %mieæ mA, mno¿ê wszystko przez 1e3, czyli mno¿e ostatecznie przez 1e-1
    figure(5)
%        semilogx( Na/1e6,log10(exp(Ea_Io_v3(2,:))) )
        loglog( Na/1e6,exp(Ea_Io_v3(2,:))*q*1e-1 )
%        semilogx( Na/1e6,log10(exp(Ea_Io(2,:))) )
%        title('log10(Io) ~ domieszkowania')
    title('Io[mA/cm2] ~ domieszkowania')
    

    tmp=find(-Ea_Io_v3(1,:)*k/q<0.095,1)';    
    figure(6)
        semilogy(1./Temperature', p_L2(:,tmp).*all_mobilities(:,tmp).*exp(-q*Vb(:,tmp)./...
        (k*Temperature')).*Emax(Vb(:,tmp),0,Na(tmp),epsR*eps),'x','MarkerSize',12)
        hold on
        y_est = polyval(Ea_Io_v3(:,tmp),1./Temperature);
        semilogy(1./Temperature,exp(y_est)','--','LineWidth',2)
        hold off
        title('Przykladowy Arrhenius')
%     x = polyfit(1./Temperature', log( p_L2(:,44).*uc.*exp(-q*Vb(:,44)./(k*Temperature')) ),1)

    % Ea z doœwiadczenia to oko³o 0.875
    fprintf('\nPoni¿ej wyszukuje pierwszy Ea, który jest mniejszy ni¿ 0.095\n');
    fprintf('Indeks elementu, ktorego Ea jest najbli¿ej 0.085eV:\n %d \n',tmp);
    fprintf('Energia aktywacji w poboli¿u tego elementu: \n %f  %f %f %f %f\n',(-Ea_Io_v3(1,(tmp-2:tmp+2))*k/q));
    fprintf('pr¹d nasycenia (Io/Sq) w pobli¿u dla tego elementu: \n %e %e %e %e %e\n',exp(Ea_Io_v3(2,(tmp-2:tmp+2))));
    fprintf('domieszkowanie dla tego elementu [cm-3]: \n %e \n',Na(tmp)/1e6);

%(-c(1)*k/q)
%exp(c(2))
end
%% FUNKCJE
   function ni = ni_fun(T,e_eff_mass,h_eff_mass,Eg,k)
        %efektywna gêstoœæ stanów pas przewodnictwa w cm-3
        Nc = 2.5*10^19*((e_eff_mass)^(3/2))*(T/300)^(3/2);
        %efektywna gêstoœæ stanów pas walencyjny w cm-3
        Nv= 2.5*10^19*((h_eff_mass)^(3/2))*(T/300)^(3/2);
        % gêstoœæ noœników w przewodnik samoistnym w m-3
        ni = sqrt(Nv*Nc)*exp(-Eg/(2*k*T))*1e6;
   end
   function Vb = Vb_fun(Na,W)
        q    = 1.6021766208*10^(-19);
        eps  = 8.85*10^-12;%A^2 s^4 / kh m^3 
        epsR = 13;
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
   function [pgb,p_L2,W,Vb,Ef] = sove_problem_for_all_Na(Na_x,pgb_x,p_L2_x,W_x,Vb_x,Ef_x,Fermi_graniczny_x,ni,T,delta,L,Qt,Et,delta_Et,Eg,epsR,deltaE,q,uc,ugb)
        i=1;
        is_fermi_small_enough = true;
        fprintf('Lewa strona\n')
        while i <= length(Na_x) && is_fermi_small_enough == true
            %fprintf('%f\n', i)
            Ef_x(i)   = findFermi(Na_x(i), T, delta, L, Qt, Et, delta_Et, Eg, epsR);
            p_L2_x(i) = p_L2_fun(ni,Ef_x(i),T);
            W_x(i)    = L/2;
            Vb_x(i)   = Vb_fun(Na_x(i),W_x(i));
            pgb_x(i)  = pgb_fun(p_L2_x(i),Vb_x(i),deltaE,T);

            if (Vb_x(i)>0.13 && Vb_x(i)<0.15) || (Vb_x(i)>0.08 && Vb_x(i)<0.1)
                fprintf('Vbi:   %f    u:   %e    Na:   %e\n', Vb_x(i),u_fun(uc, ...
                                                              Fgb_fun(delta,L,uc,ugb,p_L2_x(i),pgb_x(i)), ...
                                                              Fc_fun(W_x(i),L,delta,Vb_x(i),T))*1e4, Na_x(i)/1e6)
            end

            if Ef_x(i) <= Fermi_graniczny_x(i) || -Ef_x(i) > 0.5*Eg/q 
                is_fermi_small_enough = false;
            end
            i=i+1;

        end
        fprintf('N*    %e\n',Na_x(i-2)/1e6)
        if is_fermi_small_enough == false
            fprintf('Prawa strona\n')
            i=i-1;
            is_fermi_small_enough = true;
            while i <= length(Na_x) && is_fermi_small_enough == true
                Ef_x(i)   = Fermi_graniczny_x(i);
                p_L2_x(i) = Na_x(i);%p_L2_fun(ni,Ef_x(i),T);
                W_x(i)    = findW(Na_x(i), T, delta, L, Qt, Et, delta_Et, Eg,epsR);
                Vb_x(i)   = Vb_fun(Na_x(i),W_x(i));
                pgb_x(i)  = pgb_fun(p_L2_x(i),Vb_x(i),deltaE,T);

                if (Vb_x(i)>0.13 && Vb_x(i)<0.15) || (Vb_x(i)>0.08 && Vb_x(i)<0.1)
                fprintf('Vbi:   %f    u:   %e    Na:   %e\n', Vb_x(i),u_fun(uc, ...
                                                              Fgb_fun(delta,L,uc,ugb,p_L2_x(i),pgb_x(i)), ...
                                                              Fc_fun(W_x(i),L,delta,Vb_x(i),T))*1e4, Na_x(i)/1e6)
                end

                if -Ef_x(i) > 0.5*Eg/q 
                    is_fermi_small_enough = false;
                end
                i=i+1;

            end
        end

       pgb = pgb_x;
       p_L2 = p_L2_x;
       W = W_x;
       Vb = Vb_x;
       Ef = Ef_x;
   end
   function [Ea, all_mobilities] = calculate_activation_energy(max_Na,T,W,L,delta,Vb,uc,ugb,p_L2,pgb)
    
    all_mobilities = zeros(length(T),max_Na);
    for i=1:length(T)
       all_mobilities(i,1:max_Na)= u_fun(uc, ...
                                         Fgb_fun(delta,L,uc,ugb,p_L2(i,1:max_Na),pgb(i,1:max_Na)), ...
                                         Fc_fun(W(i,1:max_Na),L,delta,Vb(i,1:max_Na),T(i)));
    end
    Ea = zeros(2,max_Na);
    for i=1:max_Na
        Ea(:,i) = polyfit(1./T,log( all_mobilities(:,i)' ),1);
    end
     q=1.6021766208*10^(-19); %C
     k = 1.38064852*10^(-23); %J/K
     Ea(1,:) = -Ea(1,:)*k/q;
   end
%     function W = W_fun(eps,epsR,Vb,Na)
%         q=1.6021766208*10^(-19); W = sqrt(2*eps*epsR*Vb.*(q*Na).^-1);
%     end   
%    function Q_pgb = Q_pgb_fun(delta,p_L2,Vb,deltaE,T)
%         Q_pgb = delta * pgb_fun(p_L2,Vb,deltaE,T);
%    end  
%    function Q_dep_reg = Q_dep_reg_fun(Na,L)
%         Q_dep_reg = Na*L;
%    end      
%    function part_of_Qt = part_of_Qt_fun(alpha,delta_eps)
%         k = 1.38064852*10^(-23);
%         T = 300;
%         mnoznik = (k*T/(sqrt(pi)*delta_eps));
%         part_of_Qt =  mnoznik.*(1.64696 - 0.844494 * exp((0.386513 - 0.0743295 .* alpha).* alpha) + ... 
%           exp((-0.386513 - 0.0743295 .*alpha) .*alpha) .* (0.844494 - 0.844494 * erf(0.416573 - 0.46392.*alpha)) + ...
%           0.844494 * exp((0.386513 - 0.0743295 .*alpha) .*alpha).* erf(0.416573 + 0.46392 .* alpha) - ... 
%           1.64696 * erf(0.5381 .* alpha));
%    end   
%    function Q_trap = Q_trap_fun(Qt,alpha1,alpha2,delta_eps)
%         Q_trap    = Qt*(part_of_Qt_fun(alpha1,delta_eps)-part_of_Qt_fun(alpha2,delta_eps));
%    end 
  
%     function uncompensated_charge = charge_neut_cond(Ef,Vb,Na,L,Qt,delta,T,alpha1,alpha2,delta_eps)
%          k = 1.38064852*10^(-23);
%          q=1.6021766208*10^(-19);
%          e_eff_mass = 0.09;
%          h_eff_mass = 0.72;
%          Nc = 2.5*10^19*((e_eff_mass)^(3/2))*(T/300)^(3/2);
%          Nv= 2.5*10^19*((h_eff_mass)^(3/2))*(T/300)^(3/2);
%          Eg = 1.2*q;
%          ni = sqrt(Nv*Nc)*exp(-Eg/(2*k*T))*1e6;
%          deltaE = 0.05*q;
%          uncompensated_charge = Na*L -  ...
%                                 Qt*(part_of_Qt_fun(alpha1,delta_eps)-part_of_Qt_fun(alpha2,delta_eps))-...
%                                 delta*ni*exp(-q*Ef/(k*T))*exp(-(q*Vb-deltaE)/(k*T));
%     end
%% Kod przydatny, acz ju¿ nie wykorzystywany

    %tmp = find(-Ea_Io_v3(1,1:(length(Ea_Io_v3)-20))*k/q>0.085)';
    %tmp=max(tmp);

% 
% fileID = fopen('mobility_100K.txt','w');
% data = [Na/1e6;u_fun(uc, ...
%                 Fgb_fun(delta,L,uc,ugb,p_L2,pgb), ...
%                 Fc_fun(W,L,delta,Vb,T))*1e4];
% fprintf(fileID,'%f %f\n',data);
% fclose(fileID);


% loglog(Na/1e6,u_fun(uc, ...
%                 Fgb_fun(delta,L,uc,ugb,p_L2,pgb), ...
%                 Fc_fun(W,L,delta,Vb,T))*1e4)


%     fprintf('Gêstoœæ noœników na granicy ziarna:  %e cm^-3\n', pgb_fun(p_L2,Vb,deltaE,T)/1e6);
%     fprintf('Warstwa zubo¿ona:  %f um\n', W_fun(eps,epsR,Vb,Na)*1e6);
%     fprintf('Fc: %f\n', Fc_fun(W,L,delta,Vb,T));
%     fprintf('Fgb: %f\n', Fgb_fun(delta,L,uc,ugb,p_L2,pgb));
%     fprintf('u = %f cm^2/Vs\n', u_fun(uc,Fgb,Fc)*1e4);
    
%     alpha1  = (q*Ef(i)+Et+q*Vb(i))/(k*T);
%     alpha2  = (Et)/(k*T);
%     fprintf('charge:    %e,    charge_warzub   %e,   charge_trap   %e,   charge_gb   %e\n' ...
%             , charge_neut_cond(Ef(i),Vb(i),Na(i),L,Qt,delta,T,alpha1,alpha2,delta_eps) ...
%             , Q_dep_reg_fun(Na(i),L)   ...
%             , Q_trap_fun(Qt,alpha1,alpha2,delta_eps) ...
%             , Q_pgb_fun(delta,p_L2(i),Vb(i),deltaE,T) );
%     fprintf('Q2   %f, Q1   %f\n',part_of_Qt_fun(alpha2,delta_eps), part_of_Qt_fun(alpha1,delta_eps));


% testy
% delta_ef = 0.015*q;
% alpha1  = (q*Ef(1)+delta_ef+Et+q*Vb_fun(Na(1),L/2))/(k*T);
% alpha2  = (Et)/(k*T);    
% fprintf('charge:    %e,     dQt:   %e,    Q1:   %f,    Q2:   %f\n', ...
%          charge_neut_cond(Ef(1)+delta_ef,Vb_fun(Na(1),L/2),Na(1), ...
%                           L,Qt,delta,T,alpha1,alpha2,delta_eps), ...
%           Q_trap_fun(Qt,alpha1,alpha2,delta_eps), ...
%           part_of_Qt_fun(alpha1,delta_eps), ...
%           part_of_Qt_fun(alpha2,delta_eps));     
% charge_neut_cond(0,Vb(1),Na(1),L,Qt,delta,T,alpha1,alpha2,delta_eps)