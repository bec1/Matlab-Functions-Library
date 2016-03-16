function [ outp ] = imagedata_tiff( folderpath, varargin )
    
%% Extract list of fits images
listing = dir(fullfile(folderpath,'*.fits'));
imnames = cellfun( @(x) x(1:end-5), {listing(:).name}', 'UniformOutput', false);
impaths = cellfun( @(x) fullfile(folderpath,x), {listing(:).name}', 'UniformOutput', false);

new_impaths = cellfun( @(x) fullfile(folderpath,'converted',[x,'.tiff']), imnames, 'UniformOutput', false );
mkdir(fullfile(folderpath,'converted'));

%% Convert all images and save the results
for i = 1:length(new_impaths)
    temp = imagedata(imnames{i});
    imwrite(temp.od, new_impaths{i});
end

end
