function [ k_U, P_U, U_i, EF_U, n_U, index1, index2 ] = EoS_nz_KvsP_s( z_i, n_z, varargin )
%% Information
% Extract from n(z) the kappa(U) and P(U)
%
% Inputs are z_i and n_z
%   z_i : z_i == 0 refers to the trap center
%       : must be in SI units
%   n_z : density in SI units
%
% Outputs are k_U, P_U, U_i
%   U_i : the binned potential values
%   P_U : P/P0 in SI units
%   k_U : kappa/kappa0 in SI units
%
% Name value pairs include bins, trap omega, pixel, plot
% 

%% Constants
% Universal Constants
uconst.h = 6.62607004e-34;
uconst.hbar = uconst.h / (2*pi);
uconst.massLi6 = 9.988346e-27;

% Experimental Constants, CHANGE ACCORDINGLY WITH THE EXPERIMENT
econst.trapw = 2*pi*23.9;

% Other variables
bindz = 0.5;
methods = {'sum','diff',0, 2}; % {sum or trapz or simps, diff or poly, length and order for poly}
plotset = {1,1, 3, 50}; % {?main with k and P, ?technical details plot, U max (kHz), z max (um) }
Urange = [0.1 3];

% Process inputs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'trap omega', econst.trapw = varargin{i+1};
        case 'plot', plotset = varargin{i+1};
        case 'Urange', Urange = varargin{i+1};
        case 'bins', bindz = varargin{i+1};
        case 'methods', methods = varargin{i+1};
    end
end

%% Procedure

% Create potential U(z)
U_z = 0.5*uconst.massLi6*econst.trapw^2*z_i.^2;

% Sort n_z and U_z in increasing U_z with variable names U_i and n_U
[sU_z, sIndex] = sort(U_z);
sn_z = n_z(sIndex);
sz_i = z_i(sIndex);

% Dynamic Binning
U_i = zeros(size(sU_z));
n_U = zeros(size(sU_z));
z_U = zeros(size(sU_z));
dz = z_i(2) - z_i(1);
bindz = 0.5;

i=1; binIndex = 1;
while i<=length(sU_z)
    members = [i];
    j = i+1;
    while j<=length(sU_z) && (abs(sz_i(j))-abs(sz_i(i))) < (dz*bindz)
        members = [members, j]; j = j+1;
    end
    U_i(binIndex) = mean(sU_z(members));
    n_U(binIndex) = mean(sn_z(members));
    z_U(binIndex) = mean(abs(sz_i(members)));
    i = j; binIndex = binIndex + 1;
end
U_i = U_i(1:binIndex-1);
n_U = n_U(1:binIndex-1);
z_U = z_U(1:binIndex-1);

% Calculate fermi energy
EF_U = real(uconst.hbar^2 / (2*uconst.massLi6) * (6*pi^2*n_U).^(2/3));

% Calculate P/P0 and k/k0.
P_U = zeros(size(EF_U));
k_U = zeros(size(EF_U));

% Calculate P_U
switch methods{1}
    case 'sum', for i = 2:length(P)-1, for j = i:length(U_i)-1, P(i) = P(i) + (U_i(j+1)-U_i(j))*n_U(j); end; end
    case 'trapz', for i = 2:length(P)-1, P(i) = trapz(U_i(i:end),n_U(i:end)); end
    case 'simps', for i = 2:length(P)-1, P(i) = simps(U_i(i:end),n_U(i:end)); end
end
P = P ./ (2/5*n_U.*EF_U);
P = P(2:end-1); 

% Calculate k_U
switch methods{2}
    case 'diff', k_U(1:end-1) = - diff(EF_U) / (U_i(2)-U_i(1));
    case 'poly' 
        warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
        for i = 2:length(k_U)-1
            t = methods{3};
            p = polyfit(U_i(max(1,i-t):min(length(U_i),i+t)), EF_U(max(1,i-t):min(length(U_i),i+t)), methods{4});
            k_U(i) = - polyval(polyder(p),U_i(i));
        end
end

% Limit the range of U for P_U and k_U
if ~isempty(Urange)
    Urange = Urange * uconst.h * 1e3;
    index1 = find(U_i > Urange(1),1);
    index2 = find(U_i < Urange(2)); index2 = index2(end);
else
    index1 = 1;
    index2 = length(P_U);
end

%% Figure
if plotset{1}
    figure;
    subplot(2,2,1);
    grid on; title('Density and Potential vs z');  xlim([z_i(1),z_i(end)]*1e6); xlabel('z (\mum)');
    yyaxis left
    plot(z_i*1e6,n_z,'.');  ylabel('n (m^{-3})');
    yyaxis right
    plot(z_i*1e6,U_z/(uconst.h*1e3),'.'); ylabel('U (kHz)');

    subplot(2,2,2);
    grid on;  title('Density and Fermi energy vs U'); xlabel('U (kHz)'); xlim([0,plotset{3}]); 
    yyaxis left
    plot(U_i/(uconst.h*1e3),n_U,'.');  ylabel('n (m^{-3})');
    yyaxis right
    plot(U_i/(uconst.h*1e3),EF_U/(uconst.h*1e3),'.');  ylabel('E_F (kHz)');

    subplot(2,2,3);
    grid on;  title('Compressibility and Pressure vs U');  xlabel('U (kHz)'); xlim([0,plotset{3}]);
    yyaxis left
    plot(U_i(1:end-1)/(uconst.h*1e3),k_U,'.'); ylim([0 5]);  ylabel('\kappa / \kappa_0');
    yyaxis right
    plot(U_i(1:end-1)/(uconst.h*1e3),P_U,'.'); ylim([0 5]); ylabel('P / P_0');

    subplot(2,2,4);
    plot(P_U(index1:index2),k_U(index1:index2),'r.'); xlim([0 4]); ylim([0 4]);
    title('\kappa / \kappa_0 vs P / P_0'); xlabel('P / P_0'); ylabel('\kappa / \kappa_0'); grid on;
end

if plotset{2}
    figure;
    
    subplot(2,2,1);
    grid on; title('Dynamic Binning'); xlim([0 plotset{4}]); xlabel('z [\mum]');
    yyaxis left
    plot()
    plot(sU_z/(uconst.h*1e3), sn_z, 'r.', bU/(uconst.h*1e3), bn, 'bo'); xlim([0 plot_maxU]);
    
end

end

