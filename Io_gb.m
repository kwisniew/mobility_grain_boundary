function Io_gb = Io_gb(p_L2,mi_eff,Vb,Temperature,Na,epsR,eps)
    q=1.6021766208*10^(-19);
    k = 1.38064852*10^(-23); %J/K
    Io_gb = (1e-1*q*p_L2.*mi_eff.*exp(-q*Vb./...
              (k*Temperature')).*Emax(Vb,0,Na,epsR*eps))';
end