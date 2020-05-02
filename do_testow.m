fis = 0.07;
A = 2e-6;
c = polyfit(1./Temperature', ...
            log(A*uc*Nv(Temperature).*...
            Emax(fis,0,1e21,epsR*eps).*...
            exp(-q*fis./(k*Temperature)) )'...
            ,1);
figure(6)
plot(1./Temperature', log(A*uc*Nv(Temperature).*...
                      Emax(fis,0,1e21,epsR*eps).*...
                      exp(-q*fis./(k*Temperature)) )')
(-c(1)*k/q)
exp(c(2))

c = polyfit(1./Temperature',log(all_mobilities(:,tmp)),1);
plot(1./Temperature', log(all_mobilities(:,tmp)))

%% Do zapamiêtania
    %ZMIENNE
    %temperatura
    % UWAGA: ka¿da zmiana T powoduje, ¿e
    %        trzeba jeszcze raz przeliczyæ ca³ki w Mathematice!!!
    Temperature = [300,250,200,150,100];
    T=Temperature(5);%K
    %szerokoœæ ziarna
    L=1e-6;%m
    %szerokoœæ granicy ziarna
    delta=2e-9;%m
    % gêstoœæ przestrzenna ³adunku na granicy ziaren 
    Qt=0.2e10*1e4;%m^-2
    % po³o¿enie Et wzglêdem Ei
    Et=0.25*q;%eV
    % szerokoœæ po³ówkowa 
    % UWAGA: ka¿da zmiana delta_Et (bazowo: 0.083q) powoduje, ¿e
    %        trzeba jeszcze raz przeliczyæ ca³ki w Mathematice!!!
    delta_Et = 0.083*q;
    %przerwa energetyczna
    Eg = 1.20*q;
    % przenikalnoœæ elektryczna wzglêdna
    epsR = 13;