%% ImagingResonance
%  Find the imaging resonance from a set of images
%


%% Syntax
%  ImagingResonance(images)
%  ImagingResonance(images,varargin,crop)
%
%  [nums,freqs,clouds,imgresfit] = ImagingResonance(images,crop)

%% Description
%  ImagingResonance(images) produces a figure displaying the resonance from
%  a sequence of images provided in the cell array of strings input.
%  ImagingResonance(images,crop) finds the resonance from images and crops
%  each image using crop. 

%% Examples

% Read an image set

ImagingResonance(images);


