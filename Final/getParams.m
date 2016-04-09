function params = getParams(images,paramname)
%% GETPARAMS gets parameters from a sequence of images
%           Usage:   params = getParams(top_images, 'RF23')
%           Inputs:     images: a cell of strings for filenames
%                       paramname: a string of the parameter name
%           Output:  params: a cell array of parameters for each image

%% Format filenames
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
    if isempty(strfind(fullpath,'fits'))
        endpoint = length(fullpath);
    else 
        endpoint = length(fullpath)-5;
    end
    filenames{i} = fullpath(startpoint:endpoint);
end

%% Get the parameters from the snippet files
for i=1:length(filenames)
        img_name = filenames{i};
        snipout = GetSnippetValues(img_name,{paramname},'SnippetFolder','R:\Snippet');
        params{i} = str2double(snipout.value{1});
end

end