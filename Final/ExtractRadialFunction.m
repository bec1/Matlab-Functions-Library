function [ radfun, rads ] = ExtractRadialFunction( dat, center, rads, dr )

%% check for nargin
if nargin == 1
    center = Gauss2DFit(dat);
    rads = 1:length(dat(:,1))/4;
    dr = 1;
elseif nargin == 2
    rads = 1:length(dat(:,1))/4;
    dr = 1;
elseif nargin == 3
    dr = rads(2) - rads(1);
elseif nargin > 4
    disp('WARNING! TOO MANY INPUTS');
end

%% Calculate the integral
inttot = RadialIntegral_temp(dat,center,rads(1)-0.5*dr);
radfun = rads;
for i = 1:length(rads)
    integ = RadialIntegral_temp(dat,center,rads(i)+0.5*dr);
    radfun(i) = integ - inttot;
    inttot = integ;
end
radfun = radfun ./ (2*pi*dr*rads);
end

