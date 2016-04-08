function fouts = findLocalImages(filenames,varargin)
%% FINDLOCALIMAGES finds images locally
%       Usage:  fouts = findLocalImages(filenames,sourcepath)
%       Inputs:   - filenames: a cell array of string filenames           
%                 - sourcepath: a string of the source containing 2016 etc
%       Outputs:  - fouts: a cell array of full paths to local files

%% Specify the source path
switch nargin
    case 1
       sourcepath = ('\\Elder-pc\j\Elder Backup Raw Images');
    case 2
        sourcepath = varargin{1};
end

%% Find the images
    for i=1:length(filenames)
        filename = filenames{i};
        year = filename(7:10);
        month = filename(1:2);
        day = filename(4:5);
        if str2num(filename(12:13))<7 %if we're doing super late night stuff, we've slipped back in time
            day = num2str(str2num(day)-1,'%02.f');
        end
        yrmnth = strcat(year,'-',month);
        yrmnthday = strcat(yrmnth,'-',day);
        filepath = strcat(sourcepath,'\',year,'\',yrmnth,'\',yrmnthday);
        fouts{i} = strcat(filepath,'\',filename,'.fits');
    end
end