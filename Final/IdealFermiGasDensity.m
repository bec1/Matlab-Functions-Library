function [ outp ] = IdealFermiGasDensity( y_vec_micron, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program creates density profiles for the ideal Fermi gas at finite
% temperature along the axial direction of the hybrid trap. For checking 
% purposesVarious EOS quantities are calculated as function of the axial
% direction.
% In addition gaussian white noise can be added to the simluated data.
% Input parmeters are listed direclty below.


%% Input parameters

% Select T/TF and density in center of trap
TTilde = 0.1;
n0 = 10^(15); % in 1/m^3

% Select axial trapping frequency
omega_y = 2*pi * 23.9;

% y_vec_micron
% y_vec_micron = -200:0.1:200;

% make plots
makeplot = [0 0];

% Signal 2 Noise
s2n = 50;

% physical constants
hbar = 1.0545718*10^(-34);
kB = 1.38064852*10^(-23);

% Lithium parameters
mLi = 9.9883414*10^(-27);

%% Process inputs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'T', TTilde = varargin{i+1};
        case 'plot', makeplot = varargin{i+1};
        case 'n0', n0 = varargin{i+1};
        case 's2n', s2n = varargin{i+1};
    end
end

%% Generate EOS data for a large regime
Xstart = -5;
Xstop = 5;
X_vec = linspace(Xstart,Xstop, 5000);
Z_vec = exp(X_vec);

TTildeAll = 4*pi./(6*pi^2*(-PolyLogFrac(3/2,-Z_vec)).^(2/3));

%% Find fugacity for given T/TF in center of trap
[diff,SelectIndex] = min(abs(TTildeAll - TTilde));
Z_select = Z_vec(SelectIndex);

%% Find absolute temperature for given density in the center of the trap
TAbsolut = (n0/(-PolyLogFrac(3/2,-Z_select)))^(2/3) * (2*pi*hbar^2/(kB*mLi));

%% Chemical potential in the center of the trap
mu0=log(Z_select)*kB*TAbsolut;

%% chemical potential in local density approximation for hybrid box in pixel steps
y_vec = y_vec_micron*10^(-6);
mu_simulated = mu0 - 1/2 * mLi * omega_y^2 * (y_vec).^2;

%% simulated fugacity
Z_simulated = exp(mu_simulated/(kB*TAbsolut));

%% simulated density distribution
n_simulated = (mLi*kB*TAbsolut/(2*pi*hbar^2))^(3/2)*(-PolyLogFrac(3/2,-Z_simulated));

if makeplot(1)
    figure(1)
    plot(y_vec_micron,n_simulated)
    xlabel('\mu m')
    ylabel('density [1/m^3]')
end

%% PTilde as a function of y_vec
PTilde = 10*pi/(6*pi^2)^(2/3) * ...
    (-PolyLogFrac(5/2,-Z_simulated)./(-PolyLogFrac(3/2,-Z_simulated)).^(5/3));

%% KappaTilde as a function y_vec 
KappaTilde = (6*pi^2)^(2/3)/(6*pi) * ...
    (-PolyLogFrac(1/2,-Z_simulated)./(-PolyLogFrac(3/2,-Z_simulated)).^(1/3));

%% add gaussian noise on density distribution
n_simulated_noise = n_simulated + max(n_simulated) / s2n * randn(size(n_simulated));

if makeplot(2)
    figure(2)
    plot(y_vec_micron,n_simulated_noise)
    xlabel('\mu m')
    ylabel('density [1/m^3]')

    save('polarized_simulated_T0_1_noise_hires.mat','n_simulated','y_vec',...
        'n_simulated_noise','PTilde','KappaTilde')
end

%% Gather outputs
outp.z = y_vec;
outp.n = n_simulated;
outp.n2 = n_simulated_noise;
outp.Pt = PTilde; 
outp.kt = KappaTilde;


end

