function diffout = centraldiff_2(pointsinput,datainput)
%% CENTRALDIFF_2 calculates central finite differences of data(points)
%           First derivative with accuracy 2
%           Inputs:
%               data: the data to be differentiated
%               points: the x grid. Can be of nonuniform grid spacing
% 

%% Input handling

if(length(datainput) ~= length(pointsinput))
    msgbox('Please enter data of equal lengths')
end

data = reshape(datainput,[1 length(datainput)]);
points = reshape(pointsinput, [1 length(datainput)]);


%% Calculate diffs

diffout = data;

for i=1:length(data)
    if i==1 %forward difference at edge
        spacingfactor = 1/(points(2)-points(1));
        diffcoeffs = spacingfactor.*[-1, 1];
        neighborhood = data(i:i+1);
    elseif i==length(data) %backward difference at edge
        spacingfactor = 1/(points(end)-points(end-1));
        diffcoeffs = spacingfactor.*[-1 1];
        neighborhood = data(i-1:i);
    else %central difference at center
        dx = diff(points);
        spacingfactor = [-dx(i)/(dx(i-1)*(dx(i)+dx(i-1))), (dx(i)-dx(i-1))/(dx(i)*dx(i-1)),  dx(i-1)/(dx(i)*(dx(i)+dx(i-1)))];
        diffcoeffs = spacingfactor;
        neighborhood = data(i-1:i+1);
    end
    
    diffout(i) = diffcoeffs*neighborhood';
end

%% Output handling
% make output the same shape as the input
diffout = reshape(diffout,size(pointsinput));
        