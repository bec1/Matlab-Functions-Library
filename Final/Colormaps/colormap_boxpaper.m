function cmap = colormap_boxpaper(points)
    % If points is not provided, then set it to 100
    if nargin < 1, points = 100; end
    
    % Load a function for RGB's (generated from 11 points from http://colorbrewer2.org)
    RBG = load('colormap_boxpaper.mat');
    
    % Generate colormap matrix of points x 3
    points_list = linspace(0,1,points); points_list = points_list';
    cmap = [RBG.fitRs(points_list), RBG.fitBs(points_list), RBG.fitGs(points_list)];
    
end