function Emax = Emax(fi,V,Na,es)
    q=1.6021766208*10^(-19);
    Emax = sqrt(2*q.*(fi-V).*Na/es);
end