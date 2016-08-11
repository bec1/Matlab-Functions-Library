function [nums,clouds] = AtomNumberFluctuations(images)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMAGINGRESONANCE finds the resonance from a sequence of images
%
% Usage:
%
%   - Input: images: a cell of filenames, without the extension
%            crop: the usual crop parameters passed to imcrop [x0 y0 w l]
%
%   - Output: ImagingResonance(images)
%       nums = # of atoms, a.u.
%       freqs = imaging frequencies extracted from the snippet Server
%       clouds = images loaded
%
%%

images = processPaths(images);
clouds = getClouds(images);
nums = getNums(clouds);

plot(fliplr(nums),'.','MarkerSize',15)
ylim([0 1.1*max(nums)])
% hold all
% ax = plot(imgresfit);
% set(ax,'DisplayName',strcat('\nu_0 = ',num2str(imgresfit.x0)))
% xlabel('Frequency [MHz]')
% ylabel('# [a.u.]')
set(gca,'FontSize',14)
grid on

end


function imgout = processPaths(images)
    c=clock;
    year = num2str(c(1));
    month = strcat(year,'-',num2str(c(2),'%02d'));
    day = strcat(month,'-',num2str(c(3),'%02d'));
    sourcepath = '\\Elder-pc\j\Elder Backup Raw Images';
    source=strcat(sourcepath,'\',year,'\',month,'\',day);
    imgout = cellfun(@(x) strcat(source,'\',x,'.fits'), images,'UniformOutput',false);
end

function nums = getNums(clouds)
    for i=1:size(clouds,3)
        nums(i) = sum(sum(clouds(:,:,i)));
    end
    
end


function clouds = getClouds(filenames)
%% crop depending on side or top image
    if isempty(strfind(filenames{1},'top'))
        cropper = {'rect',111 ,420,150,150};
    else
        cropper = {'rect',275 ,243,150,400};
    end

%% Get the data
 data = imagedata_list(filenames,'crop',cropper);
 
 for i=1:length(data)
     clouds(:,:,i) = data(i).od2;
 end
 
 
end


function rf = getFreqs(images)

addpath('C:\Users\BEC1\Documents\GitHub\Data-Explorer-GUI\Snippet-Functions')
%% Get filenames
for i=1:length(images)
    fullpath = images{i};
    f1 = findstr(fullpath,'/');
    f2 = findstr(fullpath,'\');
    f=[f1 f2];
    if isempty(f)
        startpoint = 1;
    else
        startpoint = 1+max(f);
    end
    filenames{i} = fullpath(startpoint:end-5);
end

for i=1:length(filenames)
        img_name = filenames{i};
        snipout = GetSnippetValues(img_name,{'Side green evap'},'SnippetFolder','R:\Snippet');
        rf{i} = str2double(snipout.value{1});
end
end


