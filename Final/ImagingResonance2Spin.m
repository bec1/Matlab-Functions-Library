function [nums1,freqs1,nums2,freqs2] = ImagingResonance2Spin(images,varargin)
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

% Separate the images into the two spin states
[images1, images2] = separateImages(images);

% % Get the proper paths
% images1 = processPaths(images1)'
% images2 = processPaths(images2)'

% Read the clouds
clouds1 = getClouds(images1);
clouds2 = getClouds(images2);

% Get the atom numbers
nums1 = getNums(clouds1);
nums2 = getNums(clouds2);

% Get the imaging frequencies
freqs1 = cell2mat(getFreqs(images1));
freqs2 = cell2mat(getFreqs(images2));


%% Get the resonances
figure;
imgresfit1 = imagingResFit(freqs1,nums1);
plot(freqs1,nums1,'.','MarkerSize',15)
xlim([min(freqs1)-1 max(freqs1)+1])
hold all
ax = plot(imgresfit1);
set(ax,'DisplayName',strcat('\nu_0 = ',num2str(imgresfit1.x0)))
xlabel('Frequency [MHz]')
ylabel('# [a.u.]')
set(gca,'FontSize',14)
grid on

figure;
imgresfit2 = imagingResFit(freqs2,nums2);
plot(freqs2,nums2,'.','MarkerSize',15)
xlim([min(freqs2)-1 max(freqs2)+1])
hold all
ax = plot(imgresfit2);
set(ax,'DisplayName',strcat('\nu_0 = ',num2str(imgresfit2.x0)))
xlabel('Frequency [MHz]')
ylabel('# [a.u.]')
set(gca,'FontSize',14)
grid on


end

function [images1, images2] = separateImages(images)

    k=1; l=1;
    for i=1:length(images)
        if ~isempty(strfind(images{i},'TopA'))
            images1(k) = images(i);
            k = k+1;
        elseif ~isempty(strfind(images{i},'TopB'))
            images2(l) = images(i);
            l=l+1;
        end
    end

end

function imgresfit = imagingResFit(freqs,nums)
    guess_freq = mean(freqs);
    [xData, yData] = prepareCurveData( freqs, nums );
    ft = fittype( 'a+ b/((x-x0)^2+c)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.435698684103899 10000 0.923379642103244 guess_freq];

    % Fit model to data.
    [imgresfit, ~] =  fit( xData, yData, ft, opts );

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
    if ~isempty(strfind(filenames{1},'A'))||~isempty(strfind(filenames{1},'B'))
        cropper = {'rect',1054,851,200,200};
    elseif isempty(strfind(filenames{1},'top'))
        cropper = {'rect',110 ,110,200,200};
    else
        cropper = {'rect',272 ,182,150,300};
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
        if ~isempty(strfind(images{i},'TopA'))
            paramname = 'ImagFreq1';
        elseif ~isempty(strfind(images{i},'TopB'))
            paramname = 'ImagFreq2';
        end
        snipout = GetSnippetValues(img_name,{paramname},'SnippetFolder','R:\Snippet');
        rf{i} = str2double(snipout.value{1});
end
end


