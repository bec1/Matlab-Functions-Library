function plotsequence(filenames,varargin)
%% PLOTSEQUENCE plots a sequence of images with a varied parameter shown in the x axis
%   Usage:    plotsequence(filenames,'Wait Time')
%                   where filenames = cell array of filenames
%                   and   'Wait Time' = string of List-bound variable


clouds = getClouds(filenames);

%% Get the parameter varied and sort
s = size(clouds);
switch nargin
    case 1
        params = 1:s(3);
    case 2
        params = cell2mat(getParams(filenames,varargin{1}));
end

%% Get the images
[params,ix] = sort(params);
clouds = clouds(:,:,ix);

%% Plot the sequence
figure; imagesc(reshape(clouds,size(clouds(:,:,1),1),[],1));

set(gca,'XTick',s(2)*(1:length(params))-s(2)/2,'XTickLabel',num2str(params,'%2.2f\n'),'FontSize',14,'YTick',[]);

axis image

end

function clouds = getClouds(filenames)
%% crop depending on side or top image
    if isempty(strfind(filenames{1},'top'))
        cropper = {'rect',110 ,110,200,200};
    else
        cropper = {'rect',262 ,242,30,60};
    end

%% Get the data
 data = imagedata_list(filenames,'crop',cropper);
 
 for i=1:length(data)
     clouds(:,:,i) = data(i).od2;
 end
 
 
end
