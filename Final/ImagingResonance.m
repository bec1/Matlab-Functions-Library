function [nums,freqs,clouds,imgresfit] = ImagingResonance(images,varargin)
%% Usage:
%   - Input: images: a cell of filenames, without the extension
%            crop: the usual crop parameters passed to imcrop [x0 y0 w l]
%   - Output: ImagingResonance(images)
%       nums = # of atoms, a.u.
%       freqs = imaging frequencies extracted from the snippet Server
%       clouds = images loaded
%%
switch nargin
    case 1
        crop=[20,234,150,150];  
    case 2
        crop = varargin{1};
end
        

images = processPaths(images);

data = loadDataset(images);
clouds = getClouds(data,crop);

nums = getNums(clouds);
freqs = cell2mat(getFreqs(images));

if range(freqs)==0 
    disp('Please vary the imaging frequency')
    return
end

imgresfit = imagingResFit(freqs,nums);
plot(freqs,nums,'.','MarkerSize',15)
xlim([min(freqs) max(freqs)])
hold all
ax = plot(imgresfit);
set(ax,'DisplayName',strcat('\nu_0 = ',num2str(imgresfit.x0)))
xlabel('Frequency [MHz]')
ylabel('# [a.u.]')
set(gca,'FontSize',14)
grid on

end

function imgresfit = imagingResFit(freqs,nums)
    [xData, yData] = prepareCurveData( freqs, nums );
    ft = fittype( 'a+ b/((x-x0)^2+c)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.435698684103899 10000 0.923379642103244 180];

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

function clouds = getClouds(data,crop)
%% Get 1D axial profiles
    for i=1:length(data)
        clouds(:,:,i) = imcrop(data(i).img,crop);
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
        snipout = GetSnippetValues(img_name,{'ImagFreq1'});
        rf{i} = str2double(snipout.value{1});
end
end

function data = loadDataset(img_list)
%% LOAD_DATASET loads the raw images as OD arrays
    % Initialize data struct
    data(1:length(img_list)) = struct('name','','img',[]);
    % Load the images from the filenames
    fprintf('\n');
    for i =1:length(img_list)
        if ~isempty(img_list{i})
            fprintf('.');
            data(i).name = img_list{i};
            data(i).img=loadFitsimage(data(i).name);
        end
    end
    fprintf('\n');
end

function img = loadFitsimage(filename)
 data=fitsread(filename);
    absimg=(data(:,:,2)-data(:,:,3))./(data(:,:,1)-data(:,:,3));

%     % Identify "burned pixels" and make sure they will come out zero.
%     burnedpoints = absimg < 0;
%     absimg(burnedpoints) = 1;
% 
%     % Same thing for points which should be accidentally equal to zero
%     % (withatoms == background) or infinity (withoutatoms == background)
%     zeropoints = absimg == 0;
%     absimg(zeropoints) = 1;
% 
%     infpoints = abs(absimg) == Inf;
%     absimg(infpoints) = 1;
% 
%     nanpoints = isnan(absimg);
%     absimg(nanpoints) = 1;

%replace the pixels with a value of negtive number,0 or inf or nan by the
%average of nearset site.
    ny=size(absimg,1);
    nx=size(absimg,2);
    burnedpoints = absimg <= 0;
    infpoints = abs(absimg) == Inf;
    nanpoints = isnan(absimg);
    Change=or(or(burnedpoints,infpoints),nanpoints);
    NChange=not(Change);
    for i=2:(ny-1)
        for j=2:(nx-1)
            if Change(i,j)
                n=0;
                rp=0;
                if NChange(i-1,j)
                    rp=rp+absimg(i-1,j);
                    n=n+1;
                end
                if NChange(i+1,j)
                    rp=rp+absimg(i+1,j);
                    n=n+1;
                end
                if NChange(i,j-1)
                    rp=rp+absimg(i,j-1);
                    n=n+1;
                end
                if NChange(i,j+1)
                    rp=rp+absimg(i,j+1);
                    n=n+1;
                end
                if (n>0)
                    absimg(i,j)=(rp/n);
                    Change(i,j)=0;
                end
            end
        end
    end
    absimg(Change)=1;
    img = log(absimg);
end