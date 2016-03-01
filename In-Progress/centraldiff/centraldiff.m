function diffout = centraldiff(points,data,varargin)
%% FDIFF calculates central finite differences of data(points)
%           Inputs:
%               data: the data to be differentiated
%               points: the x grid
%               order: the order of differentiation
%               accuracy: order of qtruncation errors

%% Input handling

if(length(data) ~= length(points))
    msgbox('Please enter data of equal lengths')
end

data = reshape(data,[1 length(data)]);
points = reshape(points, [1 length(data)]);

switch nargin
    case 2
        order = 1;
        accuracy = 2;
    case 3
        order = varargin{1};
        accuracy = 2;
    case 4
        order = varargin{1};
        accuracy = varargin{2};
end

%% Calculate diffs

gridspacing = points(2) - points(1);

switch accuracy
    case 2
        diffcoeffs = [-1/2, 0, 1/2];
    case 4
        diffcoeffs = [1/12, -2/3, 0, 2/3,-1/12];
end
paddeddata = padarray(data,[0 3],'replicate');
diffout = data;

for i=1:length(data)
    center = i+3;
    switch accuracy 
        case 2
            neighborhood = paddeddata(center-1:center+1);
        case 4
            neighborhood = paddeddata(center-2:center+2);
    end
    diffout(i) = diffcoeffs*neighborhood' / gridspacing;
end
        