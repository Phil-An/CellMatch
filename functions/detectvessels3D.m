function [Vfiltered, Histostack_EQ]= detectvessels3D(Histostack, maxvesselsize)
% This function detects blood vessels in histology images. Blood vessels 
% must be dark tubular structures ins the histology volume.
%
% Inputs:
%           Histostack : a 3D matrix of a laser confocal microscopy 
%                        histology scan of brain tissue. 
%
%           maxvesselsize : parameter for detecting the largest diameter of
%                           blood vessels. This parameter depends on the 
%                           size of blood vessels in the images, and hence 
%                           on the type of objective used for acquisition.
%
% Outputs:
%           Vfiltered : 3D matrix of detected blood vessels.
%
%           Histostack_EQ : a histogram-equalized and pre-processed version
%                           of the histology laser confocal microscopy scans
%
%
% Function is written by Philip Anner (2020)

origsize = size(Histostack);
Histostack = imresize(Histostack,0.5);

if ~isa(Histostack,'double')
    Histostack=im2double(Histostack);
end

Histostack_EQ = zeros(size(Histostack));
for i=1:size(Histostack,3)
   Histostack_EQ(:,:,i)= imadjust(Histostack(:,:,i));
   Histostack_EQ(:,:,i)= adapthisteq(Histostack(:,:,i));
end

for j=1:size(Histostack_EQ,3)
    Histostack_EQ(:,:,j) = imgaussfilt(Histostack_EQ(:,:,j),2);
end

% morphological filtering for blood vessels
maxstack = double(zeros(size(Histostack_EQ)));
strength = 2:2:maxvesselsize;
strength = [strength, 1 3 5 7 9];
for k=1:length(strength)
        deg = [0 45 90 125 145 190 225 270 300 325];
    for j=1:length(deg)
        se = strel('line',strength(k), deg(j));
        tempimg = imbothat(Histostack_EQ,se);
        maxstack = max(tempimg,maxstack);
    end
end
options = struct('FrangiScaleRange', [1 8], 'FrangiScaleRatio', 2, 'FrangiAlpha', 0, 'FrangiBeta', 1000, 'FrangiC', 50, 'verbose',true,'BlackWhite',false);
[Vfiltered,~,~,~,~]=FrangiFilter3D(maxstack,options);

Vfiltered = imresize(Vfiltered,[ origsize(1) origsize(2)]);
Vfiltered=mat2gray(Vfiltered);
Histostack_EQ = imresize(Histostack_EQ,[ origsize(1) origsize(2)]);

clear maxstack tempimg strength deg Vx Vy Vz
end


