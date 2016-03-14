function [ outp ] = imagedata_list( fnames, varargin )
%% Information
% This function takes in a list of imagenames and outputs the data in a
% structrued array. Optionally, user can provide the fields to include in
% the output and settings for imagedata function.

%% Define Variabels and Process Inputs



% Defaults for the input
sample_image = 0; % 0 or 1, image is picked randomly from the list
output_fields = {'all'}; % Fields to output
Nsat = Inf;
crop_set = {'none',10,10,4,4}; % none or rect or ellipse and [centX, centY, width or radius1, height or radius2]
bg_set = {'none',10}; % none or avg or linear, width
plot_set = {0};
pol_set = {'none',0.05}; % none or sub -- What percent of light is the wrong polarization

% Inputs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'crop', crop_set = varargin{i+1};
        case 'plot', plot_set = varargin{i+1};
        case 'Nsat', Nsat = varargin{i+1};
        case 'bg', bg_set = varargin{i+1};
        case 'pol', pol_set = varargin{i+1};
        case 'sample', sample_image = varargin{i+1};
        case 'fields', output_fields = varargin{i+1};
    end
end


%% Procedure

% Prepare
outp = struct([]);
total_images = length(fnames);

% Sample
if sample_image
    temp = imagedata(fnames{randi(total_images,1)},'crop',crop_set,'plot',{1},'Nsat',Nsat,'bg',bg_set,'pol',pol_set);
end

% Import the list of images
for i = 1:total_images
    temp = imagedata(fnames{i},'crop',crop_set,'plot',plot_set,'Nsat',Nsat,'bg',bg_set,'pol',pol_set);
    if strcmp(output_fields{1},'all')
        
        for fn = fieldnames(temp)'
            outp(i).(fn{1}) = temp.(fn{1});
        end
    else
        for j = 1:length(output_fields)
            outp(i).(output_fields{j}) = temp.(output_fields{j});
        end
    end
end


%% Gather Outputs



end

