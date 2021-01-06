function [cropped_stack] = crop_stack(inputstack, rect)
% This function crops a 3D imagestack crop_stack(STACK, RECT)
%   Input: 
%           inputstack : 3D matrix to be cropped
%
%           rect : rectangle [x, y, width, height] defines croppping 
%                  coordinates 
% Function is written by Philip Anner (2020)
for i=1:size(inputstack,3)
    cropped_stack(:,:,i) = imcrop(inputstack(:,:,i), rect);
end
end


