function Out = FiniteD( x,y,SD )
%FINITED Summary of this function goes here
%   Detailed explanation goes here
N=length(x);
Out=y*0;

for i=1:N
    k1=max(1,i-SD);
    k2=min(N,i+SD);
    Out(i)=(y(k2)-y(k1))/(x(k2)-x(k1));
end

end

