function diffout = polydiff(points,data,varargin)
%% POLYDIFF does a LSQ polyfit in [-range,range] and evaluates the derivative
% Usage: polydiff(points,data,range,degree)

%% Handle inputs
switch nargin
    case 2
        range = 5;
        degree = 2;
    case 3
        range = varargin{1};
        degree = 2;
    case 4 
        range = varargin{1};
        degree = varargin{2};
end

if(length(data) ~= length(points))
    msgbox('Please enter data of equal lengths')
end

data = reshape(data,[1 length(data)]);
points = reshape(points, [1 length(data)]);

%% Evaluate polynomial derivative
N = length(data);

for i=1:N
    k1=max(1,i-range);
    k2=min(N,i+range);
    neighborx=points(k1:k2);
    neighbory=data(k1:k2);
    p=polyfit(neighborx,neighbory,degree);
    q=polyder(p);
    diffout(i)=polyval(q,points(i));
end