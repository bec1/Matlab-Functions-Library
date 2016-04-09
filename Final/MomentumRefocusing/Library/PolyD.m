function dy = PolyD( x,y,order )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p=polyfit(x,y,order);
dp=polyder(p);
dy=polyval(dp,x);

end

