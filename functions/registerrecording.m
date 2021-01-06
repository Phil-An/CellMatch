function [REC] = registerrecording(REC, varargs)
% This function co-registers all images from an in vivo ca2+ recording and
% returns a structure with aligned images.
% Pre-processing operations are defined by Multireg (VARARGS 1)
%
%
% Inputs:
%           REC : struct containing all elements of an in vivo ca2+
%           recording
%
%           Multireg (VARARG 1) : Struct for pre-processing prior to image
%           registration procedure. Multireg.preprocess may contain
%           following boolean variables:
%           -preprocess.crop : crop all in vivo acquired images defined by
%                       rect prior to all other preprocessing operations
%                       (Multireg.preprocess.crop = 1
%                       Multireg.preprocess.rect = [x, y, width, height])
%                       registerrecording(REC, Multireg)
%           -preprocess.rot90 : rotate all images clockwise for 90°
%                       (Multireg.preprocess.rot90 = 1)
%           -preprocess.fliplr : flip all images
%                       (Multireg.preprocess.fliplr= 1)
%           -preprocess.matchimagesize : automatically increase (pad) the
%                       size of all images to the same size as the image in
%                       VARARGS(2) - (RefImage)
%                       (Multireg.preprocess.matchimagesize= 1)
%                       registerrecording(REC, Multireg, RefImage)
%           -postprocess.crop : crop all in vivo acquired images defined by
%                       rect after to all other preprocessing operations
%                       (rotation, cropping, ...)
%                       (Multireg.postprocess.crop = 1
%                       Multireg.postprocess.rect = [x, y, width, height])
%                       registerrecording(REC, Multireg, RefImage)
%           -preprocess.register : co-register all in vivo acquired images
%                       using transformation parameters stored in
%                       REC.O_trans and REC.Spacing
% Function is written by Philip Anner (2020)

%% preprocess Ca2+ images
preprocess.rot90 = 0;
preprocess.fliplr = 0;
preprocess.crop = 0;
preprocess.matchimagesize = 0;
preprocess.postcropCa = 0;
preprocess.register = 0;

% check if transformation parameters are defined
if isfield(REC,'O_trans')
    preprocess.register = 1;
end

% Evaluate activated pre-processing operations
Multireg=[];
if exist('varargs')
    Multireg = varargs(1);
    %rotation
    if isfield(Multireg.preprocess, 'rot90') && Multireg.preprocess.rot90 >0
        for i=1:Multireg.preprocess.rot90
            % individual cell images
            preprocess.rot90 = preprocess.rot90 + 1;
        end
    end
    
    %flipping
    if isfield(Multireg.preprocess, 'fliplr') && Multireg.preprocess.fliplr > 0
        preprocess.fliplr = 1;
    end
    
    
    %cropping
    if isfield(Multireg.preprocess, 'cropca')
        if Multireg.preprocess.cropca == 1
            preprocess.crop = 1;
        end
    end
    
    % match image size
    if isfield(Multireg.preprocess, 'matchimagesize')
        if Multireg.preprocess.matchimagesize == 1
            preprocess.matchimagesize = 1;
            t = varargs(2);
            preprocess.Historefimage = t(:,:,1);
        end
    end
    
    % check for post processing
    if isfield(Multireg,'postprocess')
        if isfield(Multireg.postprocess, 'crop')
            preprocess.postcropCa = 1;
        end
    end
end

% Perform pre-processing operations for all elements of a recording
if isfield(REC,'IC')
    REC.IC = preprocess_element(REC.IC,preprocess, Multireg)
    if preprocess.register
        for i=1:size(REC.IC,3)
            IC(:,:,i)=bspline_transform(REC.O_trans,REC.IC(:,:,i),REC.Spacing, 3);
        end
        REC.IC=IC;
    end
end

if isfield(REC,'areas')
    REC.areas = preprocess_element(REC.areas,preprocess, Multireg)
    if preprocess.register
        for i=1:size(REC.areas,3)
            areas(:,:,i)=bspline_transform(REC.O_trans,REC.areas(:,:,i),REC.Spacing, 3);
        end
        REC.areas=areas;
    end
end


if isfield(REC,'soma')
    REC.soma = preprocess_element(REC.soma,preprocess, Multireg)
    if preprocess.register
        for i=1:size(REC.soma,3)
            soma(:,:,i)=bspline_transform(REC.O_trans,REC.soma(:,:,i),REC.Spacing, 3);
            t=regionprops(soma(:,:,i),'Centroid');
            if isempty(t)
                center=[0 0];
            else
                center=deal(t(:).Centroid);
            end
            REC.centroids(:,:,i) = center;
            perim(:,:,i)=bwperim(soma(:,:,i));
        end
        REC.soma = soma;
        REC.perim = perim;
    end
end

if isfield(REC,'vesselimg')
    REC.vesselimg = preprocess_element(REC.vesselimg,preprocess, Multireg)
    if preprocess.register
        REC.vesselimg = bspline_transform(REC.O_trans, REC.vesselimg, REC.Spacing, 3);
    end
end

if isfield(REC,'vesselimg_bw')
    REC.vesselimg_bw = preprocess_element(REC.vesselimg_bw,preprocess, Multireg)
    if preprocess.register
        REC.vesselimg_bw = bspline_transform(REC.O_trans, REC.vesselimg_bw, REC.Spacing, 3);
    end
end


if isfield(REC,'DSA')
    REC.DSA = preprocess_element(REC.DSA,preprocess, Multireg)
    if preprocess.register
        REC.DSA = bspline_transform(REC.O_trans, REC.DSA, REC.Spacing, 3);
    end
end

if isfield(REC,'minrec')
    REC.minrec = preprocess_element(REC.minrec,preprocess, Multireg)
    if preprocess.register
        REC.minrec = bspline_transform(REC.O_trans,REC.minrec,REC.Spacing, 3);
    end
end

if isfield(REC,'allIC')
    REC.allIC = preprocess_element(REC.allIC,preprocess, Multireg)
    if preprocess.register
        REC.allIC = bspline_transform(REC.O_trans,REC.allIC,REC.Spacing, 3);
    end
end

if isfield(REC,'maxrec')
    REC.maxrec = preprocess_element(REC.maxrec,preprocess, Multireg)
    if preprocess.register
        REC.maxrec = bspline_transform(REC.O_trans,REC.maxrec,REC.Spacing, 3);
    end
end

if isfield(REC,'meanrec')
    REC.meanrec = preprocess_element(REC.meanrec,preprocess, Multireg)
    if preprocess.register
        REC.meanrec = bspline_transform(REC.O_trans,REC.meanrec,REC.Spacing, 3);
    end
end

if isfield(REC,'std')
    REC.std = preprocess_element(REC.std,preprocess, Multireg)
    if preprocess.register
        REC.std = bspline_transform(REC.O_trans,REC.std,REC.Spacing, 3);
    end
end
if isfield(REC,'stdrec')
    REC.stdrec = preprocess_element(REC.stdrec,preprocess, Multireg)
    if preprocess.register
        REC.stdrec = bspline_transform(REC.O_trans,REC.stdrec,REC.Spacing, 3);
    end
end

if isfield(REC,'dff') && isfield(REC.dff,'minrec')
    REC.dff.minrec = preprocess_element(REC.dff.minrec,preprocess, Multireg)
    if preprocess.register
        REC.dff.minrec = bspline_transform(REC.O_trans,REC.dff.minrec,REC.Spacing, 3);
    end
end

if isfield(REC,'dff') && isfield(REC.dff,'maxrec')
    REC.dff.maxrec = preprocess_element(REC.dff.maxrec,preprocess, Multireg)
    if preprocess.register
        REC.dff.maxrec = bspline_transform(REC.O_trans,REC.dff.maxrec,REC.Spacing, 3);
    end
end
if isfield(REC,'dff') && isfield(REC.dff,'meanrec')
    REC.dff.meanrec = preprocess_element(REC.dff.meanrec,preprocess, Multireg)
    if preprocess.register
        REC.dff.meanrec = bspline_transform(REC.O_trans,REC.dff.meanrec,REC.Spacing, 3);
    end
end
if isfield(REC,'dff') && isfield(REC.dff,'stdrec')
    REC.dff.stdrec = preprocess_element(REC.dff.stdrec,preprocess, Multireg)
    if preprocess.register
        REC.dff.stdrec = bspline_transform(REC.O_trans,REC.dff.stdrec,REC.Spacing, 3);
    end
end
end


function Image = preprocess_element(Image,preprocess, Multireg)

if preprocess.crop
    Image = imcrop3D(Image,Multireg.preprocess.rect);
end

if preprocess.rot90 > 1
    for i=1:preprocess.rot90
        for imgindex = 1:size(Image,3)
            Image(:,:,imgindex) = rot90(Image(:,:,imgindex));
        end
    end
end

if preprocess.fliplr
    for imgindex = 1:size(Image,3)
        Image(:,:,imgindex) = fliplr(Image(:,:,imgindex));
    end
end

if preprocess.matchimagesize
    for imgindex = 1:size(Image,3)
        Image(:,:,imgindex) = matchimagesize(Image(:,:,imageindex),preprocess.Historefimage);
    end
end

if preprocess.postcropCa
    Image = imcrop3D(Image,Multireg.postprocess.rect);
end

end

