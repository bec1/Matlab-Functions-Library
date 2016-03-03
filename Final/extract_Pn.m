function [ n_i, P_n ] = extract_Pn( nxy, varargin )
%% Information
% Computes P(n) dn = # of atoms with density between n and n+dn
% 
% Inputs 
%   nxy : density of atoms in SI (1/m^3)
%
% Outputs
%   n_i : binned density vector in SI (1/m^3)
%   P_n : # of atoms with density n
% 
% Optional Name Value pairs
%   bins, nmin, nmax, plot

%% Constants and Parameters
% Universal Constants
uconst.h = 6.62607004e-34;
uconst.hbar = uconst.h / (2*pi);
uconst.massLi6 = 9.988346e-27;

% Experimental Constants, CHANGE ACCORDINGLY WITH THE EXPERIMENT
econst.trapw = 2*pi*23.9;
econst.px = 10e-6;

% Other variables
bins = 50;
nmin = min(nxy(:));
nmax = max(nxy(:));
create_plot = 1;
plot_axis = nan;
plot_title = 'P(n)';
plot_nxy = 1;
plot_nxy_pos = [.7 .7 .2 .2];
n_i = [1];

% Process inputs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'bins', bins = varargin{i+1};
        case 'nmin', nmin = varargin{i+1};
        case 'nmax', nmax = varargin{i+1};
        case 'plot', create_plot = varargin{i+1};
        case 'plot_axis', plot_axis = varargin{i+1};
        case 'plot_title', plot_title = varargin{i+1};
        case 'plot_nxy', plot_nxy = varargin{i+1};
        case 'plot_nxy_pos', plot_nxy_pos = varargin{i+1};
        case 'n_i', n_i = varargin{i+1};
    end
end

%% Procedure
nxy_flat = nxy(:);  % Flatten the matrix

% Create bins
if n_i(1) == 1
    n_i = linspace(nmin, nmax, bins+1)'; % One additional bin for the end, will be removed later
end
P_n = zeros(length(n_i)-1,1);
n_i_edges = n_i - 0.5*(n_i(2) - n_i(1));

whichbin = discretize(nxy_flat, n_i_edges);
for i = 1:length(P_n)
    binMembers = nxy_flat(whichbin == i);
    if isempty(binMembers), binMembers = 0; end
    P_n(i) = sum(binMembers);
end
n_i = n_i(1:end-1);

% Renormalize to total atom numbers
P_n = P_n * (sum(nxy_flat)) / (sum(P_n)*(n_i(2)-n_i(1)));

%% Plots
if create_plot
    % If plot axis wasn't provided, make one.
    if ~ishandle(plot_axis)
        figure; plot_axis = subplot(1,1,1);
    end
    
    % Plot
    axes(plot_axis);
    plot(plot_axis,n_i,P_n,'b-','MarkerSize',14);
    title(plot_title);
    xlabel('n');
    ylabel('# of atoms between n and n+dn');
    grid on;
    
    if plot_nxy
        % Inset
        axes('Position',plot_nxy_pos);
        imagesc(nxy); colormap gray; colorbar;
        axis off;
        axis image;
    end
    
end


end

