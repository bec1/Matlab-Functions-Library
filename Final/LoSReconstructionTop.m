function [ outp ] = LoSReconstructionTop( filename, varargin )
%% Information
% Inputs Required
%   filename : '02-21-2016_04_04_05_top' or '02-21-2016_04_04_05_top.fits'
%               or 'C://02-21-2016_04_04_05_top.fits'
% 
% Output is a structure with following fields
%   .od_r             : reconstructed od 
%   .od               : original od cropped
%   .wa               : with atoms cropped image with dark subtracted
%   .woa              : without atoms cropped image with dark subtracted
%   .map              : 1 / line of sight integration, same size as crop image
%   .map2             : fixed map to remove Inf, i.e., map(map>10) -> 0
%   .imdata           : a structure outputed from imagedata function.
%                       fields can be found in imagedata function info.
%   .[x0,rx,ry]_fit   : fitted center radius in x and y 
%   .[x0,rx,ry]_line  : extended lines from the above fitted points
%   .xsection_area_px : cross-sectional area across the cropped image in pixel^2 units.
% 
% 
% Inputs Nave-Value pairs
%   cropset : {'type',center_x,center_y,width,height} where type could be
%               rect or ellipse.
%          default  {'rect',271,218,135,350}
%   bgset : {'method',width} where method cloud be avg or (more needs to be implemented)
%               and width usually 10 to 20 pixels
%          default  {'avg',20}
%   ycuts : range of y (as a column vector) where atoms counts are good to fit disks
%          default (120:1:215)'
%   Nsat : Saturation N counts 
%          default  Inf
%   plot : {make_plot} where make_plot = 1 or 0 (more settings need to be implemented)
%          default  {0} (no plot)
% 
% Procedure Outline
%   Load a croped od image using imagedata function
%   define a range of cuts for fitting disks
%   extened the fitted centers and radii to the entire cropped image with linear extrapolation.
%   get cross sectional area using the extrapolated radii
%   create a map such that map*od gives flat image, i.e., map = 1/line of sight integration.
% 
% 
% Last Edited
%   03/01/2016 Parth : Version 1.0 completed
% 
% Bug Report
% 
% 
% 

%% Setup variables and inputs

% Default values
cropset = {'rect',271,218,135,350};
bgset = {'avg',20};
Nsat = Inf;
ycuts = (120:1:215)';
plotset = {0};

% Inputs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'crop', cropset = varargin{i+1};
        case 'bgset', bgset = varargin{i+1};
        case 'Nsat', Nsat = varargin{i+1};
        case 'plot', plotset = varargin{i+1};
        case 'ycuts', ycuts = varargin{i+1};
    end
end

%% Procedure

% Import the image data
rawdata = imagedata(filename,'crop',cropset,'bg',bgset,'plot',{0},'Nsat',Nsat); data = rawdata.od2;
yall = 1:1:size(data,1);

% Fit circles to each cross-section and extend linearly
x0 = ones(size(ycuts));
rx = ones(size(ycuts));
ry = ones(size(ycuts));
for i = 1:length(ycuts)
    y = data(ycuts(i),:);
    x = 1:1:size(y,2);
    [x0(i),rx(i),ry(i)] = diskfit(x,y);
end
[~,~,~,fun] = diskfit(x,y);
x0l = fitLine(ycuts,x0); x0l = x0l(yall);
rxl = fitLine(ycuts,rx); rxl = rxl(yall);
ryl = fitLine(ycuts,ry); ryl = ryl(yall);

% Create Line of Sight map
lineofsightint = ones(size(data));
xsectionarea_px = pi*rxl.^2;
for i = 1:size(data,1)
    lineofsightint(i,:) = 2*fun([rxl(i),rxl(i),x0l(i)],x);
end
map = 1 ./ lineofsightint;
map2 = map; map2(abs(map2)>100) = 0;
data2 = map2.*data;

% Smoothen data2 for contour plot
data2_s = imgaussfilt(data2,2);

%% figures
if plotset{1}
    figure('Units','centimeters','Position',[5 5 30 15]);

    ax1 = subplot(3,4,[1,5,9]); 
    datarange = [0, 2*data(fix(cropset{5}/2),fix(cropset{4}/2)) ];
    imshow(data,datarange); set(ax1,'YDir','normal'); %colormap jet;
    hold on; 
    for i = 1:length(ycuts)
        plot(ax1,x0(i)+rx(i),ycuts(i),'r.','MarkerSize',12);
        plot(ax1,x0(i)-rx(i),ycuts(i),'r.','MarkerSize',12);
    end
    plot(ax1,x0l+rxl,yall,'g-');
    plot(ax1,x0l-rxl,yall,'g-');
    hold off;

    ax2 = subplot(3,4,[3,7,11]); 
    data2range = [0, 2*data2(fix(cropset{5}/2),fix(cropset{4}/2)) ];
    imshow(data2,data2range); set(ax2,'YDir','normal'); %colormap jet;

    ax3 = subplot(3,4,2);
    plot(x,data(ycuts(end),:),x,fun([rx(end),ry(end),x0(end)],x)); 
    axis off; title('Top cut');

    ax4 = subplot(3,4,6);
    plot(x,data(ycuts(end/2),:),x,fun([rx(end/2),ry(end/2),x0(end/2)],x)); 
    axis off; title('Middle cut');

    ax5 = subplot(3,4,10);
    plot(x,data(ycuts(1),:),x,fun([rx(1),ry(1),x0(1)],x)); 
    axis off; title('Bottom cut');

    ax6 = subplot(3,4,[4,8,12]);
    plot(xsectionarea_px,yall,pi*rx.^2,ycuts,'r.'); 
    title('Cross-Sectional Area'); 
    xlabel('Area'); ylabel('y (px)');
    ylim([1,length(yall)]);
end

%% outputs
outp.map = map;
outp.map2 = map2;
outp.x0_fit = x0;
outp.rx_fit = rx;
outp.ry_fit = ry;
outp.x0_line = x0l;
outp.rx_line = rxl;
outp.ry_line = ryl;
outp.od_r = data2;
outp.od = data;
outp.wa = imcrop(rawdata.wa,rawdata.croprect);
outp.woa = imcrop(rawdata.woa,rawdata.croprect);
outp.imdata = rawdata;
outp.xsection_area_px = xsectionarea_px;

end

function [x0,rx,ry, modelfun] = diskfit(x,y)
    % Fitting function ellipse
    % ((y-y0)/ry)^2 + ((x-x0)/rx)^2 = 1
    % y = ry * real(sqrt(1 - ((x-x0)/rx)^2)) + y0;
    modelfun = @(A,x) ( A(2) * real(sqrt(1 - ((x-A(3))/A(1)).^2)) );
    [x, y] = prepareCurveData( x, y );
    options = statset('Display','off','UseParallel',true);
    beta0 = [40,max(y),length(x)/2];
    beta = nlinfit(x,y,modelfun,beta0,options);
    x0 = beta(3);
    rx = beta(1);
    ry = beta(2);
end

function fl = fitLine(x,y)
    [xData, yData] = prepareCurveData( x, y );
    ft = fittype( 'poly1' );
    [fl, ~] = fit( xData, yData, ft );
end