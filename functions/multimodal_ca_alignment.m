function AlignedCa = multimodal_ca_alignment(Allcells, RECDSA , Histology, Multireg)
% This function applies transformations to in vivo acquired Ca2+ imaging 
% data.
%
% Inputs:
%           Allcells : identified cells detected in in vivo Ca2+ imaging
%                      data. For obtaining Allcells run : mergecells(H).
%
%           RECDSA : Data structure containing in vivo detected blood 
%                    vessels. Acquired by calling InVivo_DetectBloodVessels
%
%           Histology : Data structure containing histology imaging data.
%                       Blood vessels detected in histology must be 
%                       included as Histology.vessels. Run 
%                       detectvessels3D(Histostack, maxvesselsize) first.
%
%           Multireg : Data structure for pre-processing and registration
%                      control pints for multimodal image alignment between
%                      histology and in vivo acquired Ca2+ imaging data.
%
%           RECDSA : Data structure containing in vivo detected blood 
%                    vessels. Acquired by calling InVivo_DetectBloodVessels
%
% Outputs: 
%           AlignedCa : cells detected in in vivo Ca2+ imaging, aligned
%           with post hoc histology.
%
%
% Function is written by Philip Anner (2020)

%% Detect segmented blood vessels in histology 
if Multireg.viewpoint.values.rotx ~= 0 && Multireg.viewpoint.values.roty ~= 0 
    if isfield(Histology,'rot')
        Histovessels = Histology.rot.vessels;
    else
        error('Rotation optimization used but no blood vessels in rotated histology structure found! Run createrotatedstack(Histology, Multireg) first')
    end
else
    Histovessels = Histology.vessels;
end

%% preprocess Ca2+ images
%rotation
if isfield(Multireg.preprocess, 'rot90') && Multireg.preprocess.rot90 >0
    for i=1:Multireg.preprocess.rot90
        % individual cell images
        Allcells.allarea = rot90(Allcells.allarea);
        Allcells.soma = rot90(Allcells.soma);
        Allcells.image = rot90(Allcells.image);
        % include images from tailvein injection
        RECDSA.DSA = rot90(RECDSA.DSA);
        RECDSA.vesselimg = rot90(RECDSA.vesselimg);
        RECDSA.minrec= rot90(RECDSA.minrec);
        RECDSA.vesselimg_bw= rot90(RECDSA.vesselimg_bw);
    end
end

%rotation
if isfield(Multireg.preprocess, 'fliplr') && Multireg.preprocess.fliplr > 0
    Allcells.allarea = fliplr(Allcells.allarea);
    Allcells.soma = fliplr(Allcells.soma);
    Allcells.image = fliplr(Allcells.image);
    
    % include images from tailvein injection
    RECDSA.DSA = fliplr(RECDSA.DSA);
    RECDSA.vesselimg = fliplr(RECDSA.vesselimg);
    RECDSA.minrec= fliplr(RECDSA.minrec);
    RECDSA.vesselimg_bw= fliplr(RECDSA.vesselimg_bw);
end


%cropping
if isfield(Multireg.preprocess, 'cropca')
    if Multireg.preprocess.cropca == 1
        Allcells.allarea = imcrop3D(Allcells.allarea, Multireg.carec);
        Allcells.soma = imcrop3D(Allcells.soma, Multireg.carec);
        Allcells.image = imcrop3D(Allcells.image, Multireg.carec);
        
        % include images from tailvein injection
        RECDSA.DSA = imcrop(RECDSA.DSA, Multireg.carec);
        RECDSA.vesselimg = imcrop(RECDSA.vesselimg, Multireg.carec);
        RECDSA.minrec= imcrop(RECDSA.minrec, Multireg.carec);
        RECDSA.vesselimg_bw= imcrop(RECDSA.vesselimg_bw, Multireg.carec);
    end
end

% if Ca2+ imaging os bigger than histo - cop Ca2
if isfield(Multireg.preprocess, 'reduceca')
    if Multireg.preprocess.reduceca == 1
        
        D = size(Allcells.allarea(:,:,1),1) - size(Histovessels(:,:,1),1);
        Dt = round(D/2);
        Db = D-round(D/2);
        Multireg.moving(1:Dt,:)=[];
        Multireg.moving(end-Db+1:end,:)=[];
        
        Allcells.allarea(1:Dt,:,:)=[];
        Allcells.allarea(end-Db+1:end,:,:)=[];
        
        Allcells.soma(1:Dt,:,:)=[];
        Allcells.soma(end-Db+1:end,:,:)=[];
        
        Allcells.image(1:Dt,:,:)=[];
        Allcells.image(end-Db+1:end,:,:)=[];
        
        Allcells.image(1:Dt,:,:)=[];
        Allcells.image(end-Db+1:end,:,:)=[];
        
        
        
        RECDSA.DSA(1:Dt,:,:)=[];
        RECDSA.DSA(end-Db+1:end,:,:)=[];
        RECDSA.vesselimg(1:Dt,:,:)=[];
        RECDSA.vesselimg(end-Db+1:end,:,:)=[];
        RECDSA.minrec(1:Dt,:,:)=[];
        RECDSA.minrec(end-Db+1:end,:,:)=[];
        RECDSA.vesselimg_bw(1:Dt,:,:)=[];
        RECDSA.vesselimg_bw(end-Db+1:end,:,:)=[];
        RECDSA.vesselsCA(1:Dt,:,:)=[];
        RECDSA.vesselsCA(end-Db+1:end,:,:)=[];
    end
end
% match image size
if isfield(Multireg.preprocess, 'matchimagesize')
    if Multireg.preprocess.matchimagesize == 1
        
        for i=1:size(Allcells.allarea,3)
        [allarea_t(:,:,i), ~] = matchimagesize(Allcells.allarea(:,:,i),Histovessels(:,:,1));
        [soma_t(:,:,i), ~] = matchimagesize(Allcells.soma(:,:,i),Histovessels(:,:,1));
        [image_t(:,:,i), ~] = matchimagesize(Allcells.image(:,:,i),Histovessels(:,:,1));
        [mergedorigimages_t(:,:,i), ~] = matchimagesize(Allcells.image(:,:,i),Histovessels(:,:,1));
        end
        
        Allcells.allarea = allarea_t;
        Allcells.soma = soma_t;
        Allcells.image = image_t;
        Allcells.image = mergedorigimages_t;
        % include images from tailvein injection
        [RECDSA.DSA, ~ ] = matchimagesize(RECDSA.DSA,Histovessels(:,:,1));
        [RECDSA.vesselimg, ~ ] = matchimagesize(RECDSA.vesselimg, Histovessels(:,:,1));
        [RECDSA.minrec, ~ ] = matchimagesize(RECDSA.minrec , Histovessels(:,:,1));
        [RECDSA.vesselimg_bw, ~ ]= matchimagesize(RECDSA.vesselimg_bw , Histovessels(:,:,1));
    end
end

% check for post processing
if isfield(Multireg,'postprocess')
    if isfield(Multireg.postprocess, 'cropCaHisto')
        RECDSA.minrec = imcrop(RECDSA.minrec, Multireg.postprocess.cropCaHisto);
        Allcells.soma = imcrop(Allcells.soma, Multireg.postprocess.cropCaHisto);
        Allcells.image = imcrop(Allcells.image, Multireg.postprocess.cropCaHisto);
        
        % include images from tailvein injection
        RECDSA.DSA = imcrop(RECDSA.DSA, Multireg.postprocess.cropCaHisto);
        RECDSA.vesselimg = imcrop(RECDSA.vesselimg, Multireg.postprocess.cropCaHisto);
        RECDSA.minrec = imcrop(RECDSA.minrec, Multireg.postprocess.cropCaHisto);
        RECDSA.vesselimg_bw = imcrop(RECDSA.vesselimg_bw, Multireg.postprocess.cropCaHisto);
    end
end

for i=1:size(Allcells.allarea,3)
    % align soma and all area perimeters
    AlignedCa.areas(:,:,i) = bspline_transform(Multireg.O_trans_ca,Allcells.allarea(:,:,i),Multireg.Spacing_ca, 3);
    AlignedCa.soma(:,:,i) = bspline_transform(Multireg.O_trans_ca,Allcells.soma(:,:,i),Multireg.Spacing_ca, 3);
    AlignedCa.soma(:,:,i)=imbinarize(AlignedCa.soma(:,:,i));
    AlignedCa.somaperim(:,:,i)=bwperim(AlignedCa.soma(:,:,i));
    
    % calculate centroids in aligned coordinate system
    t = regionprops(AlignedCa.soma(:,:,i),'Centroid');
    AlignedCa.somacentroid(:,:,i) = [t.Centroid i];
    
    % transform original IC images
    AlignedCa.images(:,:,i) = bspline_transform(Multireg.O_trans_ca,Allcells.image(:,:,i),Multireg.Spacing_ca, 3);
end

% transform images from tailvein injection
AlignedCa.RECDSA.DSA = bspline_transform(Multireg.O_trans_ca,RECDSA.DSA,Multireg.Spacing_ca, 3);
AlignedCa.RECDSA.vesselimg = bspline_transform(Multireg.O_trans_ca,RECDSA.vesselimg,Multireg.Spacing_ca, 3);
AlignedCa.RECDSA.minrec = bspline_transform(Multireg.O_trans_ca,RECDSA.minrec,Multireg.Spacing_ca, 3);
AlignedCa.RECDSA.vesselimg_bw = bspline_transform(Multireg.O_trans_ca,RECDSA.vesselimg_bw,Multireg.Spacing_ca, 3);

% set counter for number of temporally aligned cells
AlignedCa.nmerged=Allcells.nmerged;
