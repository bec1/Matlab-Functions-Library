function [ outp ] = LoSReconstructionTop( filename, varargin )
%% Information

%% Setup variables and inputs
% Default values
cropset = {'rect',271,218,135,350};
bgset = {'avg',20};
Nsat = Inf;
y_cuts = (120:5:215)';

% Inputs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'cropset', cropset = varargin{i+1};
        case 'bgset', bgset = varargin{i+1};
        case 'Nsat', Nsat = varargim{i+1};
    end
end

%% Procedure

% Import the image data
data = imagedata(filename,'crop',cropset,'bg',bgset,'plot',{0},'Nsat',Nsat);

% Fit circles to each pixel
[xsec_area, corrected_data] = get_xsection_area(data.od2,(120:5:215)',econst.pixelsize,1);

z_i_raw = (1:size(data.od2,1))' * econst.pixelsize;

n_z = sum(data.od2,2) ./ (xsec_area / econst.pixelsize * (uconst.sigma0 / 2)) * econst.fudge;
figure; plot(z_i_raw * 10^6, n_z,'b.'); xlabel('z (um)'); ylabel('density n(z)');
fitres = createFitTFnz(z_i_raw,n_z,1);
z_i = z_i_raw - fitres.x0;
figure; plot(z_i * 10^6, n_z,'b.'); xlabel('z (um)'); ylabel('density n(z)'); grid on;
clearvars data z_i_raw corrected_data xsec_area fitres

end

