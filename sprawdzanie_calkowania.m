%Dobrze ca³kkuje jeœli Ef jest ujemne a Vb jest dodatnie
%Poni¿ej za³o¿y³em, ¿e Ef=-Eg/2, a Vb=Eg/2 +/-delta_et
k = 1.38064852*10^(-23);
q=1.6021766208*10^(-19);
T = 300;
delta_Et = 0.083*q;
b1=((-Eg/2)+(Eg/2 + 0.05*delta_Et) )/(k*T)
Q1 =  1.64696 - 0.844494 * exp((0.386513 - 0.0743295 * b1)* b1) + ... 
      exp((-0.386513 - 0.0743295 *b1) *b1) * (0.844494 - 0.844494 * erf(0.416573 - 0.46392*b1)) + ...
      0.844494 * exp((0.386513 - 0.0743295 *b1) *b1)* erf(0.416573 + 0.46392 * b1) - ... 
      1.64696 * erf(0.5381 * b1);
  
mnoznik = (k*T/(sqrt(pi)*delta_eps));

wynik = Q1*mnoznik

x = (-10*delta_Et:delta_Et/10:10*delta_Et)/(k*T);
plot(x,part_of_Qt_fun(x))

function part_of_Qt = part_of_Qt_fun(alpha)
    k = 1.38064852*10^(-23);
    q=1.6021766208*10^(-19);
    T = 300;
    delta_Et = 0.083*q;
    delta_eps = delta_Et/(2*sqrt( log(2) ));
    mnoznik = (k*T/(sqrt(pi)*delta_eps));
    part_of_Qt =  mnoznik.*(1.64696 - 0.844494 * exp((0.386513 - 0.0743295 .* alpha).* alpha) + ... 
      exp((-0.386513 - 0.0743295 .*alpha) .*alpha) .* (0.844494 - 0.844494 * erf(0.416573 - 0.46392.*alpha)) + ...
      0.844494 * exp((0.386513 - 0.0743295 .*alpha) .*alpha).* erf(0.416573 + 0.46392 .* alpha) - ... 
      1.64696 * erf(0.5381 .* alpha));
end