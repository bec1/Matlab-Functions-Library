function [nums,freqs,clouds,imgresfit] = SlicingResonance(images,freqs)
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


clouds = getClouds(images)+.08;
if iscell(freqs)
    freqs = cell2mat(freqs);
end
nums = getNums(clouds);

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
    guess_freq = mean(freqs);
    [xData, yData] = prepareCurveData( freqs, nums );
    ft = fittype( 'a+ b/((x-x0)^2+c)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [7000 -10000 0.923379642103244 guess_freq];

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
    if isempty(strfind(filenames{1},'top'))
        cropper = {'rect',190 ,190,100,100};
    else
        cropper = {'rect',262 ,252,150,300};
    end

%% Get the data
 data = imagedata_list(filenames,'crop',cropper);
 
 for i=1:length(data)
     clouds(:,:,i) = data(i).od2;
 end
 
 
end



