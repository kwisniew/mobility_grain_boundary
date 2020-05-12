function Nv = Nv(Tem)
    %na pocz¹tku przemna¿am przez 1e6, aby mieæ jednostki m^-3
    Nv = 1e6*2.5*10^19*((0.72)^(3/2)).*(Tem/300).^(3/2);
end