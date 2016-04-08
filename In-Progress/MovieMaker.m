function MovieMaker(filenames,fileout,paramname)
%% MOVIEMAKER creates a movie from a sequence of images
%          Usage:

%% Side or top cropping?
if isempty(strfind(filenames{1},'top'))
    % Side
    cropper = {'rect',110 ,110,200,200};
    fileout = strcat(fileout,'_side');
else
    % Top
    cropper = {'rect',282 ,222,250,250};
    fileout = strcat(fileout,'_top');
end

%% Get and write the movie
params = cell2mat(getParams(filenames,paramname));
data = imagedata_list(filenames,'crop',cropper);
writeToFile(data,fileout,params,paramname);

end

function writeToFile(data,fileout,params,paramname)
%% Write the movie to file
    [params,ix] = sort(params);
    data= data(ix);
    v = VideoWriter(fileout);
    v.FrameRate = 5;
    open(v);
    image = data(1).od2;
       imagesc(image);
       axis image;
       colormap parula;
       caxis([-.1 1]);
    set(gca,'nextplot','replacechildren'); 
    for k = 1:length(data) 
       image = data(k).od2;
       imagesc(image);
       axis image;
       colormap parula;
       caxis([-.1 1]);
       titlestring = strcat(paramname, ' = ', num2str(params(k),'%2.1f'));
       title(titlestring)
       frame = getframe(gcf);
       writeVideo(v,frame);
    end

    close(v);

end