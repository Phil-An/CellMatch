function [Histology,Multireg] = createrotatedstack(Histology ,Multireg )
% This function rotates histology imaging data for maximizing blood vessel
% overlap between in vivo and post hoc acquired imaging data.
%
% Inputs:
%           Histology : Data structure containing histology imaging data.
%                       histology must include 3D matrix of images of the
%                       histology imagestack, detected blood vessels
%                       Histology.vessels, acquired by
%                       detectvessels3D(Histostack, maxvesselsize), and
%                       segmented cells Histology.segmentedcells determined
%                       by annotatesegmentedcells(imagestack,Histostack_EQ)
%
%           Multireg : Data structure for pre-processing and registration
%                      control pints for multimodal image alignment between
%                      histology and in vivo acquired Ca2+ imaging data.
%                      To obatin Multireg, call findworkingdistance(Histology,RECDSA)
% Outputs:
%           Histology : Data structure containing histology imaging data.
%                       Rotation optimized histology is stored in
%                       Histology.rot
%
%           Multireg : Data structure for pre-processing and registration
%                      control pints for multimodal image alignment between
%                      histology and in vivo acquired Ca2+ imaging data.
%                      The Datastructure includes rotation parameters that
%                      optimize blood vessel overlap between in vivo and
%                      post hoc acquired images.
%
% Function is written by Philip Anner (2020)


if isfield(Multireg, 'viewpoint')
    if Multireg.viewpoint.values.rotx ~= 0 & Multireg.viewpoint.values.roty ~= 0
        % transform image stacks of all channels
        if isfield (Histology,'Ch')
            for i = 1:size(Histology.Ch,2)
                rotstack.Ch{i}.stack= im2double(Histology.Ch{i}.stack(:,:,1:Multireg.Histo.vesselimgdepth));
                rotstack.Ch{i}.stack = imwarp(rotstack.Ch{i}.stack, affine3d(makehgtform('xrotate',deg2rad(Multireg.viewpoint.values.rotx),'yrotate',deg2rad(Multireg.viewpoint.values.roty))),'Interp', 'nearest', 'FillValues', 255);
                rotstack.Ch{i}.stack(rotstack.Ch{i}.stack >1)=0;
                rotstack.Ch{i}.stack=mat2gray(rotstack.Ch{i}.stack);
            end
        elseif isfield (Histology,'ch')
            for i = 1:size(Histology.ch,2)
                rotstack.ch{i}.stack= im2double(Histology.ch{i}.stack(:,:,1:Multireg.Histo.vesselimgdepth));
                rotstack.ch{i}.stack = imwarp(rotstack.ch{i}.stack, affine3d(makehgtform('xrotate',deg2rad(Multireg.viewpoint.values.rotx),'yrotate',deg2rad(Multireg.viewpoint.values.roty))),'Interp', 'nearest', 'FillValues', 255);
                rotstack.ch{i}.stack(rotstack.ch{i}.stack >1)=0;
                rotstack.ch{i}.stack=mat2gray(rotstack.ch{i}.stack);
            end
        else
            error('Histology channels not found! Please check data structure that: Histology.Ch{1}.stack contains histology images from channel 1 etc. ')
        end 
        
        % transform blood vessels
        rotstack.vessels= im2double(Histology.vessels(:,:,1:Multireg.Histo.vesselimgdepth));
        rotstack.refvessels= imref3d(size(rotstack.vessels));
        rotstack.vessels = imwarp(rotstack.vessels, affine3d(makehgtform('xrotate',deg2rad(Multireg.viewpoint.values.rotx),'yrotate',deg2rad(Multireg.viewpoint.values.roty))),'Interp', 'nearest', 'FillValues', 255);%, 'OutputView',rotstack.refvessels);
        rotstack.vessels(rotstack.vessels >1)=0;
        rotstack.vessels=mat2gray(rotstack.vessels);
        
        % transform segmented cells
        
        rotstack.segmentedcells= im2double(Histology.segmentedcells(:,:,1:Multireg.Histo.vesselimgdepth));
        rotstack.refsegmentedcells= imref3d(size(rotstack.segmentedcells));
        rotstack.segmentedcells = imwarp(rotstack.segmentedcells, affine3d(makehgtform('xrotate',deg2rad(Multireg.viewpoint.values.rotx),'yrotate',deg2rad(Multireg.viewpoint.values.roty))),'Interp', 'linear', 'FillValues', 255);%, 'OutputView',rotstack.refsegmentedcells);
        rotstack.segmentedcells(rotstack.segmentedcells >1)=0;
        rotstack.segmentedcells(rotstack.segmentedcells<1)=0;
        rotstack.segmentedcells=imbinarize(rotstack.segmentedcells);
        for i=1:size(rotstack.segmentedcells,3)
            rotstack.segmentedcells(:,:,i)=bwareaopen(rotstack.segmentedcells(:,:,i),20);
        end
        
        %% Export rotated Histology model
        Histology.rot = rotstack;
        
        %% match image size of Multireg
        if size(Multireg.cavesselsAffinereg) ~= size(rotstack.vessels(:,:,1))
            [Multireg.cavesselsAffinereg_fit , ~ ] = matchimagesize(Multireg.cavesselsAffinereg,rotstack.vessels(:,:,1));
            Multireg.preprocess.matchimagesize = 1;
        end
        clear rotstack.refsegmentedcells rotstack.refvessels
    else
        msgbox('No rotation required!');
    end
else
    error('No rotation optimization parameters found. Please run rotation_optimization(bloodvessels, Multireg) first!')
end
end

