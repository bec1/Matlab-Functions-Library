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
h=figure();
pixellength=1.44*10^(-6);
sigma0=0.215/2*10^(-12);
Nsat=330;
ROI1 = [205,15,150,480];
sm = 2;
nbins = 100;
CropTail=1;
IfTailTailor=1;
Fudge=2.62;
D=85;
H=25;
Volume=pi/4*D^2*H*pixellength^3;
output.Volume=Volume;
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
output.Nz=n;
output.zPixel=z;
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

%% Get the density and kF from the total atom number
Ntot=sum(n);
nintrap=Ntot/Volume;
kFn=(6*pi^2*nintrap)^(1/3);
output.Ntot=Ntot;
output.nintrap=nintrap;
output.kF_num=kFn;
%% Plot n vs kz
nvskzfit = plotnvskz(kz,n);
output.n=n;
output.nvskzfit=nvskzfit;
%% Get n1d(k) vs k^2 
mu=nvskzfit.mu*1e12;
k0=sqrt(abs(mu)); %sqrt(\mu)
kz0=nvskzfit.x0*1e6;
kz=kz-kz0;
% kmin=min(kz);kmax=max(kz);
% kgrid1=linspace(kmin,kmax,51);
% [kz,n,~,~ ]=BinGrid(kz,n,kgrid1,0 );

n1dz=n/pixellength;
n1dk=(n1dz/Volume)*hbar/(m*omega);
kzsq=kz.^2;

output.kz=kz;
output.n1dofk=n1dk;
output.n1dofz=n1dz;
output.kzsq=kzsq;

%% Bin n1d(kz^2) 
% kzsqBinGrid=linspace(0,max(sqrt(kzsq)),nbins+1).^2;
kzsqBinGrid=linspace(0,max(kzsq),nbins+1);
[ kzsqBin,n1dkBin,kzsqStd,n1dkStd ] = BinGrid( kzsq,n1dk,kzsqBinGrid,0 );
% kzsqBin=kzsq;
% n1dkBin=n1dk;
kzsqBin(isnan(n1dkBin))=[];n1dkBin(isnan(n1dkBin))=[];
kzsqStd(isnan(n1dkBin))=[];n1dkStd(isnan(n1dkBin))=[];
output.kzsqBin=kzsqBin;
output.n1dofkBin=n1dkBin;
output.kzsqStd=kzsqStd;
output.n1dkStd=n1dkStd;
%% Scale kz and plot vs kz^2
subplot(2,2,3);
scatter(kzsq,n1dk,'DisplayName','Unbinned');
hold on
plot(kzsqBin,n1dkBin,'r.-','DisplayName','Binned');
line([kFn^2,kFn^2],[min(n1dk),max(n1dk)],'DisplayName','kF from N_{tot}','color','c','LineWidth',2);
hold off
xlim([0,4.5*kFn^2]);
xlabel('k_z^2 (m^{-1})');
ylabel('n_{k,1D}');
legend('show');
%% Differentiate
[fk,fkStd]=FiniteD( kzsqBin,kzsqStd,n1dkBin,n1dkStd,sm );
fk=-8*pi^2*fk;
fkStd=-8*pi^2*fkStd;
subplot(2,2,4);
kzBin=sqrt(kzsqBin);
kzFit=kzBin/k0;
[P,ffit]=FDfit(kzFit,fk);
P(3)=P(3)*k0^2;
[EF_Fit,n_Fit,kF_Fit,beta]=GetEF(P);

T=1/(beta*EF_Fit);
errorbar(kzBin,fk,fkStd);
hold on
plot(kzBin,ffit,'DisplayName',['1/\beta\mu=',num2str(1/P(2)),'T/T_F=',num2str(T)]);
line([kFn,kFn],[0,1],'DisplayName','kF from N_{tot}','color','c','LineWidth',2);
line([kF_Fit,kF_Fit],[0,1],'DisplayName','kF from Fermi-Dirac','color','r','Linewidth',2)
legend('show')
hold off
title(['Additional prefactor =',num2str(P(1))]);
xlabel('k (m^{-1})');ylabel('n(k)');
output.kzBin=kzBin;
output.Pfit=P;
output.nfit=n_Fit;
output.EF_Fit=EF_Fit;
output.kF_Fit=kF_Fit;
output.beta=beta;
output.T=T;
output.fk=fk;
output.fkStd=fkStd;
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
    xlim([min(kz) max(kz)])
    xlabel('k_z [x 10^{6} m^{-1}]')
    ylabel('n (Atom per pixel)')
    set(gca,'FontSize',14)
%     kmin=min(kz);kmax=max(kz);
%     kgrid1=linspace(kmin,kmax,51);
%     kgrid2=linspace(kmin,kmax,101);
%     kgrid3=linspace(kmin,kmax,201);
%     [k1,n1,~,~ ]=BinGrid(kz,n,kgrid1,0 );
%     plot(k1,n1,'*','DisplayName','Nbin=50');
%     [k2,n2,~,~ ]=BinGrid(kz,n,kgrid2,0 );
%     plot(k2,n2,'DisplayName','Nbin=100');
%     [k3,n3,~,~ ]=BinGrid(kz,n,kgrid3,0 );
%     plot(k3,n3,'DisplayName','Nbin=200');
    legend('show')
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


