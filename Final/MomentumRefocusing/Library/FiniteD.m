function [Out,Outerr] = FiniteD( x,xstd,y,ystd,SD )
%FINITED Summary of this function goes here
%   Detailed explanation goes here
N=length(x);
Out=y*0;
Outerr=Out;
for i=1:N
    k1=max(1,i-SD);
    k2=min(N,i+SD);
    Out(i)=(y(k2)-y(k1))/(x(k2)-x(k1));
    Outerr(i)=sqrt(((ystd(k2)-ystd(k1))/(y(k2)-y(k1)))^2+((xstd(k2)-xstd(k1))/(x(k2)-x(k1)))^2)*Out(i);
end

end

