function [P,ffit]=FDfit(k,f)
    FDfun=@(P,k) P(1)*1./(exp(P(2)*(k.^2/P(3)-1))+1);
    P0=[1,5,1];
    P=nlinfit(k,f,FDfun,P0);
    ffit=FDfun(P,k);
end