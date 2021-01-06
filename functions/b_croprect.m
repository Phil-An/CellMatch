function [ I_crop, position] = b_croprect( varargin )
% b_croprect crop image interavtively by defining a rectangle for cropping
% or define a rectangle as input parmeter. b_croprect returns the cropped
% image (I_crop) and the rectangle for cropping (position)
%
% [ I_crop position] = b_croprect(IMAGE,RECTANGLE)
%
% Function is written by Philip Anner (2020)
rot = 0;
switch nargin
    case 1
        % No cropping rectangle provided. Use default rectangle size.
        imageG = varargin{1};
        recx = 499;
        recy = 499;
        figure, imshow(imadjust(imageG));
        title('Crop image and confirm with double click')
        h = imrect(gca, [50 50 recx recy]); % create rectangle on the image
        
    case 2
        imageG = varargin{1};
        recx= varargin{2}(3);%-1;
        recy= varargin{2}(4);%-1;
        figure, imshow(imadjust(imageG));
        title('Crop image and confirm with double click')
        h = imrect(gca, [varargin{2}(1) varargin{2}(2) recx recy]);
        h.setResizable(false);

    otherwise
        error('Wring number of input arguments! Call  b_croprect(IMAGE) or  b_croprect(IMAGE,RECTANGLE)');
end

position = wait(h); % get position
position = floor(position);

I_crop = imcrop(imageG,position ); % crop image
rectangle('Position',position,'EdgeColor','r','LineWidth', 2)
delete(h)
end
