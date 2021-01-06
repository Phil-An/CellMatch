function Multireg = multimodal_alignment(Histology, Multireg,RECDSA)
% This function is used to determine the final registration parameters for
% aligning in vivo acquired Ca2+ imaging data with post hoc histology.
%
% Inputs:
%          Histovessels : a 3D matrix representing blood vessels detected
%                         in histology. Call detectvessels3D(Histostack,
%                         maxvesselsize) for detecting blood vessels.
%
%           Multireg : Multimodal registration data acquired by
%                                 findworkingdistance(Histology,RECDSA).
%
%           RECDSA : Data structure containing in vivo detected blood
%                    vessels. Acquired by calling InVivo_DetectBloodVessels
%
% Outputs:
%           Multireg : Data structure for pre-processing and registration
%                      control pints for multimodal image alignment between
%                      histology and in vivo acquired Ca2+ imaging data.
%
%
%
% Function is written by Philip Anner (2020)


if Multireg.viewpoint.values.rotx ~= 0 & Multireg.viewpoint.values.roty ~= 0
    if isfield(Histology,'rot')
        Histovessels = Histology.rot.vessels;
    else
        error('Rotation optimization used but no blood vessels in rotated histology structure found! Run createrotatedstack(Histology, Multireg) first')
    end
else
    Histovessels = Histology.vessels;
end

%% preprocess
%rotation
if isfield(Multireg.preprocess, 'fliplr') && Multireg.preprocess.fliplr > 0
    RECDSA.minrec = fliplr(RECDSA.minrec);
    RECDSA.vesselimg = fliplr(RECDSA.vesselimg);
    RECDSA.DSA = fliplr(RECDSA.DSA);
    RECDSA.vesselimg_bw = fliplr(RECDSA.vesselimg_bw);
end

%rotation
if isfield(Multireg.preprocess, 'rot90') && Multireg.preprocess.rot90 >0
    for i=1:Multireg.preprocess.rot90
        % individual cell images
        RECDSA.minrec = rot90(RECDSA.minrec);
        RECDSA.vesselimg = rot90(RECDSA.vesselimg);
        RECDSA.DSA = rot90(RECDSA.DSA);
        RECDSA.vesselimg_bw = rot90(RECDSA.vesselimg_bw);
    end
end

%cropping
if isfield(Multireg.preprocess, 'cropca')
    if Multireg.preprocess.cropca == 1
        RECDSA.minrec = imcrop(RECDSA.minrec, Multireg.carec);
        RECDSA.vesselimg = imcrop(RECDSA.vesselimg, Multireg.carec);
        RECDSA.DSA = imcrop(RECDSA.DSA, Multireg.carec);
        RECDSA.vesselimg_bw = imcrop(RECDSA.vesselimg_bw,Multireg.carec);
    end
end


% check for post processing
if isfield(Multireg,'postprocess')
    if isfield(Multireg.postprocess, 'cropCaHisto')
        RECDSA.minrec = imcrop(RECDSA.minrec, Multireg.postprocess.cropCaHisto);
        RECDSA.vesselimg = imcrop(RECDSA.vesselimg, Multireg.postprocess.cropCaHisto);
        RECDSA.DSA = imcrop(RECDSA.DSA, Multireg.postprocess.cropCaHisto);
        RECDSA.vesselimg_bw = imcrop(RECDSA.vesselimg_bw,Multireg.postprocess.cropCaHisto);
    end
end

% match image size
if isfield(Multireg.preprocess, 'matchimagesize')
    if Multireg.preprocess.matchimagesize == 1
        [RECDSA.minrec, ~] = matchimagesize(RECDSA.minrec,Histovessels(:,:,1));
        [RECDSA.vesselimg, ~] = matchimagesize(RECDSA.vesselimg,Histovessels(:,:,1));
        [RECDSA.DSA, ~] = matchimagesize(RECDSA.DSA,Histovessels(:,:,1));
        [RECDSA.vesselimg_bw,~] =matchimagesize(RECDSA.vesselimg,Histovessels(:,:,1));
    end
end

if isfield(Multireg,'postprocess')
    if isfield(Multireg.postprocess,'matchimagesize') && Multireg.postprocess.matchimagesize==1
        for i=1:size(Histovessels,3)
            [ t(:,:,i) , ~] = matchimagesize(Histovessels(:,:,i),RECDSA.minrec);
        end
        Histovessels=t;
    end
end
%% Registration
if isfield(Multireg,'viewpoint')
    if isfield(Multireg.viewpoint,'values')
        Histovesselimg= max(Histovessels(:,:,1 : end),[],3);
    end
else
    Histovesselimg= max(Histovessels(:,:,1:Multireg.Histo.vesselimgdepth),[],3);
end
if isfield(RECDSA,'vesselimg')
    Cavesselimg = RECDSA.vesselimg;
else
    error('No vesselimg in RECDSA found! Please run InVivo_DetectBloodVessels first!')
end

Multireg.fixed=Histovesselimg;
Multireg.moving=Cavesselimg;

if ~isa(Multireg.moving,'double')
    Multireg.moving=im2double(Multireg.moving);
end

if ~isa(Multireg.fixed,'double')
    Multireg.fixed=im2double(Multireg.fixed);
end

% for better visualization
Multireg.fixed=imadjust(Multireg.fixed);
Multireg.moving=imadjust(Multireg.moving);

% check if Ca2+ imaging os bigger than histo
% top to bottom
if size(Multireg.moving,1) > size(Multireg.fixed,1)
    D = size(Multireg.moving,1) - size(Multireg.fixed,1);
    Dt = round(D/2);
    Db = D-round(D/2);
    Multireg.moving(1:Dt,:)=[];
    Multireg.moving(end-Db+1:end,:)=[];
    Multireg.preprocess.reduceca = 1;
end

% if not same size - fix it
if size(Multireg.fixed,1)~= size(Multireg.moving,1) | size(Multireg.fixed,2)~= size(Multireg.moving,2)
    [Multireg.fixed, Multireg.moving] = matchimagesize(Multireg.fixed,Multireg.moving);
    Multireg.preprocess.matchimagesize = 1;
end

% init registration
Options=struct('Similarity','cc','Registration','Both','Penalty',1e-3,'MaxRef',80,'Grid',[],'Spacing',[8 8],'MaskMoving',[],'MaskStatic',[],'Verbose',0,'Points1',[],'Points2',[],'PStrength',[],'Interpolation','Linear','Scaling',[1 1]);
Multireg.Options=Options;

% check if registration points already exist
if isfield(Multireg,'movingPoints')    
    [movingPoints, fixedPoints] = cpselect(Multireg.moving,Multireg.fixed,(Multireg.movingPoints),(Multireg.fixedPoints),'Wait',true);

elseif isfield(Multireg,'preprocess')    
    [movingPoints, fixedPoints] = cpselect(Multireg.moving,Multireg.fixed,(Multireg.preprocess.movingPoints),(Multireg.preprocess.fixedPoints),'Wait',true);
else
    [movingPoints, fixedPoints] = cpselect(Multireg.moving,Multireg.fixed,'Wait',true);
end

if exist('movingPoints','var')==1
    Multireg.movingPoints=movingPoints;
    Multireg.fixedPoints=fixedPoints;
    Multireg.PStrength=ones(size(Multireg.fixedPoints,1),1);
    Multireg.PStrength=Multireg.PStrength.*0.95;
    Multireg.Options.Points1=Multireg.movingPoints;
    Multireg.Options.Points2=Multireg.fixedPoints;
    Multireg.Options.PStrength=Multireg.PStrength;
else
    sprintf('No Registration Points defined!')
end
clear movingPoints fixedPoints
Options = Multireg.Options;
Options.Points1=fliplr(Multireg.movingPoints);
Options.Points2=fliplr(Multireg.fixedPoints);

% run image registration
[Multireg.Ireg,Multireg.O_trans_ca,Multireg.Spacing_ca,~,~,~] = image_registration(Multireg.moving,Multireg.fixed ,Options);

implay(cat(3,Multireg.fixed,Multireg.Ireg))

figure, imshowpair(Multireg.fixed,Multireg.Ireg)

% Create binary Ireg image
Multireg.Ireg_bw = bspline_transform(Multireg.O_trans_ca,RECDSA.vesselimg_bw, Multireg.Spacing_ca);
end

