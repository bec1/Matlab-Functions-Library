function output = momentumfocusRV(momimages,bgimages,varargin)
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
pixellength=1.44*10^(-6);
sigma0=0.215/2*10^(-12);
Nsat=330;
ROI1 = [205,2,150,500];
sm = 2;
nbins = 100;
CropTail=1;
IfTailTailor=1;
Fudge=2.62;
D=85;
H=25;
Volume=pi/4*D^2*H*pixellength^3;

for i =1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'SM'
                sm=varargin{i+1};
            case 'Nbins' 
                nbins=varargin{i+1};
            case 'ROI1'
                ROI1=varargin{i+1};
            case 'IfTailTailor'
                IfTailTailor=varargin{i+1};
            case 'TailRange'
                Zrange=varargin{i+1};
                zmin=Zrange(1);
                zmax=Zrange(2);
                CropTail=0;
            case 'Fudge'
                Fudge=varargin{i+1};
        end
    end
end
%% Import functions and data
addpath('C:\Users\Elder\Documents\GitHub\Matlab-Functions-Library\Final');
m= 9.96e-27;
omega = 2*pi*24;
hbar = 1.05e-34;
%% Average images
Nmom=length(momimages);
momavg=0;
for i=1:Nmom
    [~,tempraw]=imagedata(momimages{i});
    Ntemp=AtomNumber(tempraw,pixellength.^2,sigma0, Nsat);
    momavg=momavg+Ntemp;
end
momavg=momavg/Nmom;

Nbg=length(bgimages);
bgavg=0;
for i=1:Nbg
    [~,tempraw]=imagedata(bgimages{i});
    Ntemp=AtomNumber(tempraw,pixellength.^2,sigma0, Nsat);
    bgavg=bgavg+Ntemp;
end
bgavg=bgavg/Nbg;

%% BG subtraction
momimg=momavg-bgavg;
momimg=momimg*Fudge;
%% Crop images
momcrop = imcrop(momimg,ROI1);


%% Get profiles
n=sum(momcrop,2)';
z=1:length(n);
n(isnan(n))=0;
if IfTailTailor
    if CropTail
        h=figure();
        scatter(z,n);
        questdlg('Now give the range for tail fitting');
        [x,y]=getpts(h);
        close(h);
        zmin=min(x);
        zmax=max(x);
    end
    n=TailTailor(n,z,zmin,zmax);
end

%% Plot the Image
figure(1);
subplot(2,2,1); imagesc(momcrop); axis image; axis off
colormap gray; caxis([0,max(momcrop(:))]);

% %% Correct for changing radius
% OUTP = LoSReconstructionTop(refimg,'cropset',{'rect',crop(1)+round(crop(3)/2),crop(2)+round(crop(4)/2),crop(3),crop(4)});
% areas = OUTP.xsection_area_px;
% output.totalatoms = sum(bgsub_profile);
% pxsize = (13e-4) / 9;
% areas = areas*(pxsize^2);
% volume = areas*pxsize;
% area_cor_profile = bgsub_profile./volume;
% n = area_cor_profile;

%% Get kz
z = pixellength *(1:length(n));
kz = m*omega*z/hbar;

output={n,kz};
%% Plot n vs kz
nvskzfit = plotnvskz(kz,n);

%% Get n1d(k) vs k^2 
mu=nvskzfit.mu*1e12;
kF=sqrt(mu);
kz0=nvskzfit.x0*1e6;
kz=kz-kz0;
n1dz=n/pixellength;
n1dk=(n1dz/Volume)*hbar/(m*omega);
kzsq=kz.^2;

%% Bin n1d(kz^2) 
kzsqBinGrid=linspace(0,max(kzsq),nbins+1);
[ kzsqBin,n1dkBin,~,~ ] = BinGrid( kzsq,n1dk,kzsqBinGrid,0 );
kzsqBin(isnan(n1dkBin))=[];n1dkBin(isnan(n1dkBin))=[];

%% Scale kz and plot vs kz^2
subplot(2,2,3);
scatter(kzsq./(kF.^2),n1dk);
hold on
plot(kzsqBin./(kF.^2),n1dkBin,'r.-');
hold off

%% Differentiate
fk=-8*pi^2*FiniteD( kzsqBin,n1dkBin,sm );
subplot(2,2,4);
kzBin=sqrt(kzsqBin);
scatter(kzBin,fk);
hold on
[P,ffit]=FDfit(kzBin,fk);
plot(kzBin,ffit,'DisplayName',['1/\beta\mu=',num2str(1/P(2))]);
legend('show')
hold off
title(['Additional prefactor =',num2str(P(1))]);
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
    kz = kz/1e6;
    plot(kz,n,'.','MarkerSize',10,'DisplayName','data')
    
    [xData, yData] = prepareCurveData( kz, n );
    % Set up fittype and options.
    ft = fittype( 'a*log(1+exp(beta*(mu-(x-x0)^2)))+d', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [30 2 0 8 4];

    % Fit model to data.
    [fitresult, ~] = fit( xData, yData, ft, opts );
    tovertf = 1/(fitresult.beta*fitresult.mu);
    hold all
    % plot
    plot(kz,fitresult(kz),'k','LineWidth',2,'DisplayName',strcat('1/(\beta\mu) = ',num2str(tovertf,'%0.2f')))
    legend('show')
    xlim([min(kz) max(kz)])
    xlabel('k_z [x 10^{6} m^{-1}]')
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

function [P,ffit]=FDfit(k,f)
    FDfun=@(P,k) P(1)*1./(exp(P(2)*((k/P(3)).^3-1))+1);
    P0=[0,5,max(k)/2];
    P=nlinfit(k,f,FDfun,P0);
    ffit=FDfun(P,k);
end