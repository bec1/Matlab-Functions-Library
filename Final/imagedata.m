function [ data, rawdata ] = imagedata( filename, varargin )
%% Information
% Inputs Required
%   filename: string of image name. It could simply be the name of image or be a complete path.
% 
% Inputs Optional via name-value pairs
%   crop : {'rect' or 'ellipse', center X, centerY, width or radiusX, height or radius Y}
%   bg   : {'avg', width for averaging}
%   plot : {0 or 1}
%   Nsat : 1234 (Saturation count for particular type of camera settings)
% 
% Outputs
%   data is a structure. The fields are listed below
%       wa, woa   : with and without atoms images with dark counts subtracted
%       abs, abs2 : absorption images wa./woa. abs2 has pixels with value
%                   (-Inf,0] and Inf fixed with average of suurouding or set to 1
%       od        : -log(abs2)
%       od2       : exists only if cropping was specified
%       
%       badpts1, badpts2, badpts3 : troble points in abs, ones that didnt
%                    get fixed by averaging neighbors and after fixing
%                    everything (should be 0)
%       bgavg     : value of the bg removed
%       imagetype : top, side_n, side_fk_3, side_fk_4
% 
% Procedure
%   Finding the image path 
%       Check if it has already been provided -> yes, then skip search
%                                             -> no, find in major then minor paths
%       If it cannot be found, throw an error.
% 
%   Load raw fits image
%       If it contains '(' then copy to temporary file, load image, and
%       delete temporary file.
% 
%   Further analysis
%       Image type -> figure out if it is a top, side (normal, fk3, fk4) or defringed image 
%       Remove dark count -> if fk4 using both dark
%                         -> if defringed, no modification is needed
%       Absorption image by simply dividing WithAtoms/WithoutAtoms 
%       Fix absorption -> since we need -log(abs) to be real, values of abs must be within 0 < abs < inf
%                      -> For points that don't satisfy this condition, try to fix them by taking an average of the neighbors.
%                      -> and if taking average doesn't fix the problem, set abs to 1 (i.e. OD -> 0)
%       OD by simply taking -log(abs2) where abs2 stands for the fixed absorption
%       Cropping -> Crop the image with provided type (none, rect, ellipse). 
%                -> If none, do not crop the image. 
%                -> Store cropped image as od2.
% 
% Last Edits
%   Parth 2/20/2016 : Version 1.0 finished
%   Parth 2/26/2016 : Version 1.1 - can read defringed files, added visualization for bad points.

%% Database of computers
% WARMING! Pc name MUST be a valid variable name. Use the following code to determine your pc name
% disp(matlab.lang.makeValidName(char(java.lang.System.getProperty('user.name'))));

% RanchoP
pcdatabase.RanchoP.master_paths = {'/Users/RanchoP/Dropbox (MIT)/BEC1/Image Data and Cicero Files/Data - Raw Images/'};
pcdatabase.RanchoP.minor_paths = {'/Users/RanchoP/Desktop';...
                                  '/Users/RanchoP/Documents';...
                                  '/Users/RanchoP/Downloads'};
pcdatabase.RanchoP.snippet_path = '/Users/RanchoP/Dropbox (MIT)/BEC1/Image Data and Cicero Files/Data - Raw Images/Snippet_output';

% Elder
pcdatabase.Elder.master_paths = {'J:\Elder Backup Raw Images';...
                                 'C:\Users\Elder\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images'};
pcdatabase.Elder.minor_paths = {'C:\Users\Elder\Desktop'};
pcdatabase.Elder.snippet_path = 'C:\Users\Elder\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images\snippet_output';

% BEC1
pcdatabase.BEC1.master_paths = {'\\Elder-pc\j\Elder Backup Raw Images';...
                                   'C:\Users\BEC1\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images'};
pcdatabase.BEC1.minor_paths = {'C:\2016-01';...
                               'C:\2016-02';...
                               'C:\Users\BEC1\Desktop'};
pcdatabase.BEC1.snippet_path = 'C:\Users\BEC1\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images\Snippet_output';

% BEC1Top
pcdatabase.BEC1Top.master_paths = {'\\Elder-pc\j\Elder Backup Raw Images';...
                                   'C:\Users\BEC1\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images'};
pcdatabase.BEC1Top.minor_paths = {'C:\2016\2016-02';...
                                  'C:\Users\BEC1\Desktop'};
pcdatabase.BEC1Top.snippet_path = 'C:\Users\BEC1\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images\Snippet_output';

% ParthPatel, home pc
pcdatabase.ParthPatel.master_paths = {'D:\Dropbox Sync\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images'};
pcdatabase.ParthPatel.minor_paths = {'C:\Users\Parth Patel\Downloads';...
                                     'C:\Users\Parth Patel\Documents';...
                                     'C:\Users\Parth Patel\Desktop'};
pcdatabase.ParthPatel.snippet_path = 'D:\Dropbox Sync\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images\Snippet_output';

% 
% biswaroopmukherjee
pcdatabase.biswaroopmukherjee.master_paths = {'\\Elder-pc\j\Elder Backup Raw Images';...
                                   'C:\Users\BEC1\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images'};
pcdatabase.biswaroopmukherjee.minor_paths = {'C:\2016-01';...
                               'C:\2016-02';...
                               'C:\Users\BEC1\Desktop'};
pcdatabase.biswaroopmukherjee.snippet_path = 'C:\Users\BEC1\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images\Snippet_output';

%% Constants, variables and inputs
% Universal constants

% Experimental constants

% Input variables
Nsat = Inf;
crop_set = {'none',10,10,4,4}; % none or rect or ellipse and [centX, centY, width or radius1, height or radius2]
bg_set = {'none',10}; % none or avg or linear, width
plot_set = {0};

% Other Variables
imagetype = 'unknown'; % 'side_n', 'side_fk_3', 'side_fk_4', 'top', 'unknown'
pcname = char(java.lang.System.getProperty('user.name')); 
pcname = matlab.lang.makeValidName(pcname);
validpc = isfield(pcdatabase,pcname); % Add line to convert name to valid variable name
filefound = 0;

% Inputs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'crop', crop_set = varargin{i+1};
        case 'plot', plot_set = varargin{i+1};
        case 'Nsat', Nsat = varargin{i+1};
        case 'bg', bg_set = varargin{i+1};
    end
end

%% Pre-Procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Common error check for type, length
validfile = ischar(filename) && length(filename) >= 6;
if ~validfile, error('Given filename failed to pass validity check.'); end

% Add .fits extension if its not there
if ~strcmp(filename(end-4:end),'.fits'), filename = [filename,'.fits']; end

% does filename contain filepath? Does Matlab know where file is?
filefound = exist(filename,'file');

% File not yet found and pc is not right
if ~filefound && ~validpc, error(['PC name ',pcname,' is not in the database. Please add it.']); end

% Find date information
imdatestr = regexp(filename,'\d\d-\d\d-\d\d\d\d','match');
imdatenum = datenum(imdatestr,'mm-dd-yyyy');
subpaths = cell(12,1);
subpaths{1} = fullfile(datestr(imdatenum,'yyyy'),datestr(imdatenum,'yyyy-mm'),datestr(imdatenum,'yyyy-mm-dd'));
subpaths{2} = fullfile(datestr(imdatenum-1,'yyyy'),datestr(imdatenum-1,'yyyy-mm'),datestr(imdatenum-1,'yyyy-mm-dd'));
subpaths{2} = fullfile(datestr(imdatenum-2,'yyyy'),datestr(imdatenum-2,'yyyy-mm'),datestr(imdatenum-2,'yyyy-mm-dd'));
subpaths{3} = fullfile(datestr(imdatenum+1,'yyyy'),datestr(imdatenum+1,'yyyy-mm'),datestr(imdatenum+1,'yyyy-mm-dd'));

% Search for image on master_paths
for i = 1:length(pcdatabase.(pcname).master_paths)
    for j = 1:length(subpaths)
        if ~filefound && exist(fullfile(pcdatabase.(pcname).master_paths{i}, subpaths{j}, filename),'file')
            filefound = 1;
            filename = fullfile(pcdatabase.(pcname).master_paths{i}, subpaths{j}, filename);
            break;
        end
    end
    if filefound, break; end
end

% Search for image on minor_paths
for i = 1:length(pcdatabase.(pcname).minor_paths)
    if ~filefound && exist(fullfile(pcdatabase.(pcname).minor_paths{i},filename),'file')
        filefound = 1;
        filename = fullfile(pcdatabase.(pcname).minor_paths{i},filename);
        break;
    end
    subfolders = dir(pcdatabase.(pcname).minor_paths{i}); subfolders = subfolders([subfolders.isdir]); subfolders = {subfolders.name};
    for j = 1:1:length(subfolders)
        if ~filefound && exist(fullfile(pcdatabase.(pcname).minor_paths{i},subfolders{j},filename),'file')
            filefound = 1;
            filename = fullfile(pcdatabase.(pcname).minor_paths{i},subfolders{j},filename);
            break;            
        end
    end
    if filefound, break; end
end

% File still not found
if ~filefound, error('The file was not found anywhere!'); end

% Store the data
data.filepath = filename;

%% Raw Image Loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If filepath contains '(' than temporary copy
if ~isempty(strfind(filename,'('))
    copyfile(filename,fileparts(userpath),'f');
    [~,name,ext] = fileparts(filename);
    temp_filename = fullfile(fileparts(userpath),[name,ext]);
    rawdata = fitsread(temp_filename);
    delete(temp_filename);
else
    rawdata = fitsread(filename);
end


%% Post-Procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine image type
[~,fname] = fileparts(filename);
if ~isempty(strfind(fname,'top')), imagetype = 'top';
elseif ~isempty(strfind(fname,'defringed')), imagetype = 'defringed';
else imagetype = 'side';
end

if strcmp(imagetype,'side')
    if size(rawdata,3) < 3, imagetype = 'unknown'; warning('This is NOT a standard absorption image. It has less than 3 layers of data. Use only the rawdata output.');
    elseif size(rawdata,3) == 3, imagetype = [imagetype, '_n'];
    elseif sum(sum(rawdata(:,:,4)))==0, imagetype = [imagetype, '_fk_3']; 
    else imagetype = [imagetype, '_fk_4']; end
end
data.imagetype = imagetype;

% remove dark count
if strcmp(imagetype,'top') || strcmp(imagetype,'side_n') || strcmp(imagetype,'side_fk_3')
    data.wa = rawdata(:,:,1) - rawdata(:,:,3);
    data.woa = rawdata(:,:,2) - rawdata(:,:,3);
elseif strcmp(imagetype,'side_fk_4')
    data.wa = rawdata(:,:,1) - rawdata(:,:,3);
    data.woa = rawdata(:,:,2) - rawdata(:,:,4);
else
    data.wa = rawdata(:,:,1);
    data.woa = rawdata(:,:,1);
end

% absorption
data.abs = data.wa ./ data.woa;
data.abs2 = data.abs;

% data.absorption should be within (0,inf), NOT (-inf,0]
badpts = data.abs <= 0 | data.abs == Inf | isnan(data.abs);
verybadpts = zeros(size(data.abs));
data.badpts = sum(badpts(:));
data.badpts2 = 0;
l1 = size(data.abs,1); l2 = size(data.abs,2);
for i = 1:l1
    for j = 1:l2
        if badpts(i,j)
            l3 = 1;
            data.abs2(i,j) = mean(mean( data.abs(max(1,i-l3):min(l1,i+l3),max(1,j-l3):min(l2,j+l3)) ,'omitnan'),'omitnan');
            while data.abs2(i,j) <= 0 || data.abs2(i,j) == Inf || isnan(data.abs2(i,j))
                l3 = l3 + 1; data.abs2(i,j) = mean(mean( data.abs(max(1,i-l3):min(l1,i+l3),max(1,j-l3):min(l2,j+l3)) ,'omitnan'),'omitnan');
                if l3 > 4, data.abs2(i,j) = 1; data.badpts2 = data.badpts2+1; verybadpts(i,j) = 1; end
            end
        end
    end
end

% od
data.od = -log(data.abs2);
data.badpts3 = sum(sum( abs(data.od)==Inf | imag(data.od)~=0  ));

% Nsat correction
if ~(Nsat == Inf || Nsat == 0), data.od = data.od + (data.woa - data.wa) / Nsat; end

% crop image
switch crop_set{1}
    case 'none'
    case 'rect'
        rect = [crop_set{2} - crop_set{4}/2, crop_set{3} - crop_set{5}/2, crop_set{4}, crop_set{5}];
        data.od2 = imcrop(data.od, rect);
        if strcmp(bg_set{1},'avg')
            rect2 = [rect(1:2) - bg_set{2}, rect(3:4) + 2*bg_set{2}];
            tempod = imcrop(data.od, rect2);
            bgtotal = sum(tempod(:)) - sum(data.od2(:));
            bgpixels = size(tempod,1) * size(tempod,2) - size(data.od2,1) * size(data.od2,2);
            data.bgavg = bgtotal / bgpixels;
            data.od2 = data.od2 - data.bgavg;
        end
    case 'ellipse'
        rect = [crop_set{2} - crop_set{4} - 2, crop_set{3} - crop_set{5} - 2, 2*crop_set{4} + 4, 2*crop_set{5} + 4];
        data.od2 = imcrop(data.od, rect);
        if strcmp(bg_set{1},'avg')
            rect2 = [rect(1:2) - bg_set{2}, rect(3:4) + 2*bg_set{2}];
            tempod = imcrop(data.od, rect2);
            bgtotal = sum(tempod(:)) - sum(data.od2(:));
            bgpixels = size(tempod,1) * size(tempod,2) - size(data.od2,1) * size(data.od2,2);
            data.bgavg = bgtotal / bgpixels;
            data.od2 = data.od2 - data.bgavg;
        end
        map = zeros(size(data.od2)); cx = size(data.od2,2)/2; cy = size(data.od2,1)/2; rx = crop_set{4}; ry = crop_set{5};
        for i = 1:size(data.od2,1)
            for j = 1:size(data.od2,2)
                if ((j-cx)/rx)^2 + ((i-cy)/ry)^2 <= 1, map(i,j) = 1; end
            end
        end
        data.od2 = data.od2 .* map;
end

%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_set{1} == 1
    range1 = [0,max(data.woa(:))];
    range2 = [0,max(max(rawdata(:,:,3)))];
    figure;
    subplot(3,4,[1,5]); imshow(data.abs2,[0,1.2]); set(gca,'YDir','normal'); title('Absorption');
    subplot(3,4,[2,6]); imshow(data.od,[0,max(data.od(:))]);set(gca,'YDir','normal'); title('OD');
    subplot(3,4,9); imshow(data.wa,range1);set(gca,'YDir','normal'); title('With Atoms');
    subplot(3,4,10); imshow(data.woa,range1);set(gca,'YDir','normal'); title('Without Atoms');
    if strcmp(imagetype,'top') || strcmp(imagetype,'side_n') || strcmp(imagetype,'side_fk_3') || strcmp(imagetype,'side_fk_4')
        subplot(3,4,11); imshow(rawdata(:,:,3),range2);set(gca,'YDir','normal'); title('Dark 1');
    end
    if strcmp(imagetype,'side_fk_4'), subplot(3,4,12); imshow(rawdata(:,:,4),range2);set(gca,'YDir','normal'); title('Dark 2'); end
    if strcmp(crop_set{1},'rect') || strcmp(crop_set{1},'ellipse')
        subplot(3,4,[3,4,7,8]); imshow(data.od2,[0,max(data.od2(:))]);set(gca,'YDir','normal'); title('Cropped OD');
        subplot(3,4,[2,6]); hold on; rectangle('Position', rect, 'EdgeColor','red'); 
        if strcmp(bg_set{1},'avg'), rectangle('Position', rect2, 'EdgeColor','blue'); end 
        hold off;    
    end
    % Drawing bad points

    [row, col] = find(badpts);
    subplot(3,4,[2,6]); hold on;
    if ~isempty(row), plot(col, row, 'g.','MarkerSize',5); end
    hold off;
    
    [row, col] = find(verybadpts);
    subplot(3,4,[2,6]); hold on;
    if ~isempty(row), plot(col, row, 'r.','MarkerSize',5); end
    hold off;
    
end

%% Collecting outputs
data.croprect = rect;

end

