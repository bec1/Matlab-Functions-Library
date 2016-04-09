function [EF,n,kF,beta]=GetEF(P)
    m= 9.96e-27;
    hbar = 1.05e-34;
    P0=P;P0(1)=1;
    FDfun=@(P,k) P(1)*1./(exp(P(2)*(k.^2/P(3)-1))+1);
    mu=P(3)*hbar^2/(2*m);
    beta=P(2)/mu;
    k0=sqrt(P(3));
    kgrid=linspace(0,10*k0,4000);
    fk=FDfun(P0,kgrid);
    Intfun=fk.*kgrid.^2/(2*pi^2);
    n=trapz(kgrid,Intfun);
    kF=(6*pi^2*n)^(1/3);
    EF=kF^2*hbar^2/(2*m);
    
end

