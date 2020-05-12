function Io_DLTS = Io_DLTS(A,uc,Temperature,fis,epsR,eps)
    q=1.6021766208*10^(-19);
    k = 1.38064852*10^(-23); %J/K
    Io_DLTS = 1e-1*q*A*uc*Nv(Temperature).*...
              Emax(fis,0,1e21,epsR*eps).*...
              exp(-q*fis./(k*Temperature));
end