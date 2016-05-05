function output = momentumfocus(momimages,bgimages,refimg,varargin)
%% Calculate momentum focused profiles
% Inputs: momimages: a cell array of T/4 momentum space focused image filenames
%         bgimages: a cell array of no atoms image filenames for background sub
%         refimg: a single string of a juicy picture in the hybrid for cone correction
%         
% Optional inputs:
%         sm: smoothing. Use 4 or 5 for good results
%         nbins: number of bins in kz/kf squared space - 100 is good
%         
% Outputs:
%         output: a struct containing the following fields:
%             kz: kz in m^{-1}
%             profile: n(kz)
%             
%             kz_sq: kz/kf squared
%             nmeans: binned densities as a function of kz_sq
%             
%             k: k/kf 
%             nofk: n(k/kf)
%             nofkfit: the fermi dirac distribution fit to the data

%%
switch nargin
    case 3
        sm = 4;
        nbins = 70;
    case 5
        sm = varargin{1};
        nbins = varargin{2};
end
%% Import functions and data
addpath('C:\Users\Elder\Documents\GitHub\Matlab-Functions-Library\Final');
m= 9.96e-27;
omega = 2*pi*24;
hbar = 1.05e-34;

%% Average images
momavg = imgAvg(momimages);
bgavg = imgAvg(bgimages);

%% Crop images
crop = [205,0,150,500];
momcrop = imcrop(momavg,crop); figure(1);subplot(2,2,1); imagesc(momcrop); axis image; axis off
colormap gray; caxis([-0.1 0.4])
bgcrop = imcrop(bgavg,crop);

%% Get profiles
raw_profile = sum(momcrop,2);
bg_profile = sum(bgcrop,2);

%% Subtract background gradients
bgfitresult = bgsmoothfit(bg_profile);
bg_fit_profile = bgfitresult(1:length(bg_profile));
bgsub_profile = raw_profile - bg_fit_profile;

%% Correct for changing radius
OUTP = LoSReconstructionTop(refimg,'cropset',{'rect',crop(1)+round(crop(3)/2),crop(2)+round(crop(4)/2),crop(3),crop(4)});
areas = OUTP.xsection_area_px;
output.totalatoms = sum(bgsub_profile);
pxsize = (13e-4) / 9;
areas = areas*(pxsize^2);
volume = areas*pxsize;
area_cor_profile = bgsub_profile./volume;
n = area_cor_profile;

%% Get kz
z = pxsize *(1:length(n));
kz = m*omega*z/hbar;

%% Plot n vs kz
nvskzfit = plotnvskz(kz,n);

%% Scale kz and plot vs kz^2
[kzsqedges,nmeans] = plotbinkz2(kz,n,nvskzfit,nbins);

%% Differentiate
[k,nofk,nofkfit] = plotnofk(kzsqedges,nmeans,sm);

%% Output results

output.kz = kz;
output.profile = n;

output.kz_sq = kzsqedges;
output.nmeans = nmeans;

output.k =k;
output.nofk = nofk;
output.nofkfit = nofkfit;
end

function [k,nofk,fitresult] = plotnofk(kzsqedges,nmeans,sm)
%% plot nofk and try to fit
    nofk = -polydiff(kzsqedges,nmeans,sm,1);
    nofk = nofk/max(nofk);%% renormalize
    k = sqrt(kzsqedges);
    figure(1);subplot(2,2,4);
    plot(k,nofk,'r.','MarkerSize',10,'DisplayName',strcat('sm=',num2str(sm)))
    xlim([0 3])
    ylim([-.2 1.1])
    hold all
    
    [xData, yData] = prepareCurveData( k(2:end), nofk(2:end) );

    % Set up fittype and options.
    ft = fittype( 'a/(1+exp(beta*(x^2-mu)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.258582251418772 18 0.593361860386109];

    % Fit model to data.
    [fitresult, ~] = fit( xData, yData, ft, opts );
    tovertf = 1/(fitresult.beta*fitresult.mu);
    plot(kzsqedges,fitresult(kzsqedges),'k','LineWidth',2,'DisplayName',strcat('T/T_F = ',num2str(tovertf,'%0.3f')) )
    legend('show')
    
    xlim([0 2])
    xlabel('k/k_F')
    ylabel('n(k)')
    set(gca,'FontSize',14)
    hold off
    
    hold off
end

function [kzedges,nmeans] = plotbinkz2(kz,n,nvskzfit,nbins)
%% Plot n vs kz^2

kz_scale = (kz-(nvskzfit.x0 *1e8))/(sqrt(nvskzfit.mu)*1e8);
kz_sq = kz_scale.^2;
figure(1);subplot(2,2,3);
plot(kz_sq,n,'.','DisplayName','raw')
hold all
[kzedges,nmeans] = binkz2(kz_sq,n,nbins);
plot(kzedges,nmeans,'k.-','LineWidth',1,'MarkerSize',10,'DisplayName','binned')
xlim([0 3])
ylabel('n [atoms/cm^{-3}/spin]')
xlabel('(k_z/k_F)^2')
set(gca,'FontSize',14)
legend('show')
hold off

% [xData, yData] = prepareCurveData( kz_sq, n );
% 
% % Set up fittype and options.
% ft = fittype( 'a*log(1+exp(beta*(mu-x)))+d', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [0.2 2 0 4 5];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );



end

function [edges,nmeans] = binkz2(kz_sq,n,nbins)
%% bin in kz^2
[~,edges,bins] = histcounts(kz_sq,nbins);
edges = edges+ (edges(2)-edges(1))/2;
for i=1:length(edges)
    binlist = [];
    for j=1:length(n)
        if bins(j) == i
            binlist = [binlist,n(j)];
        end
    end
    nmeans(i) = mean(binlist);
end
end

function fitresult = plotnvskz(kz,n)
%% Plot n vs kz
    figure(1);subplot(2,2,2);
    n = n/1e9;
    kz = kz/1e8;
    plot(kz,n,'.','MarkerSize',10,'DisplayName','data')

    [xData, yData] = prepareCurveData( kz, n );

    % Set up fittype and options.
    ft = fittype( 'a*log(1+exp(beta*(mu-(x-x0)^2)))+d', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.2 2 0 4 5];

    % Fit model to data.
    [fitresult, ~] = fit( xData, yData, ft, opts );
    tovertf = 1/(fitresult.beta*fitresult.mu);
    hold all
    % plot
    plot(kz,fitresult(kz),'k','LineWidth',2,'DisplayName',strcat('T/T_F = ',num2str(tovertf,'%0.2f')))
    legend('show')
    xlim([min(kz) max(kz)])
    xlabel('k_z [x 10^{8} m^{-1}]')
    ylabel('n [x 10^9 atoms/cm^{-3}/spin]')
    set(gca,'FontSize',14)
    hold off
end

function bgfitresult = bgsmoothfit (bg_profile)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( [], bg_profile );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 1.2338479537501e-05;

% Fit model to data.
[bgfitresult, ~] = fit( xData, yData, ft, opts );

% figure(5)
% plot(xData,yData)
% hold all
% plot(bgfitresult)
end


function odout = imgAvg(images)
    data_out = imagedata_avg(images);
    odout = data_out.raw_avg;
end