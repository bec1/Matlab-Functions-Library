function data_out = imagedata_avg(data)
%% IMGDATA_AVG averages a dataset provided in the input data

%% Handle the case where the data hasn't yet been acquired from the folder
if ~isstruct(data)
    data = imagedata_list(data);
end

%% Average the images

for i=1:length(data)
    wa_avg(:,:,i) = data(i).wa;
    woa_avg(:,:,i) = data(i).woa;
    img_avg(:,:,i) = data(i).od;
end

data_out.od_avg = mean(img_avg,3);
data_out.wa_avg = mean(wa_avg,3);
data_out.woa_avg = mean(woa_avg,3);
data_out.abs_avg = data_out.wa_avg./data_out.woa_avg;
data_out.raw_avg = -log(data_out.abs_avg);

%% Get some outputs 
data_out.filenames = {data(:).filename};
data_out.filepaths = {data(:).filepath};
data_out.imagetype = data(1).imagetype;


end