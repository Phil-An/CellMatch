% Please do not forget to download required 3rd party software described in
% README.md

% Add required folders
addpath('./')

% Load manually derived parameters for the Sample data
load('Sampledata.mat')
%% Process in vivo Ca2+ imaging data
% Align all in vivo calcium imaging sessions and track cells across
% experiments
InVivo_Align_CaImagingSessions
%% Proces Histology

% load histology images - channel 1 
imagestack = load_histology_imagestack;

% In case you have scanned histology with multiple channels, you can load 
% histology scans for channel 2 by using:
% imagestackch2 = load_histology_imagestack;
% and copying / modyfing the lines where image transformations are 
% performed on the histology image stacks.

%% Crop histology image stacks if they are of different size

% Stack 1
% check if stack 1 and 2 are of different size. If yes - define cropping
% rectangle and crop stack 1
if size(imagestack{1},1)~= size(imagestack{2},1) | size(imagestack{1},2) ~= size(imagestack{2},2)
    % Define cropping rectangle
    [~, Histology.stack{1}.rec]=b_croprect(imagestack{1}(:,:,20));
    % Use defined rectangle to crop image stack
    Histology.stack{1}.Ch{1} = crop_stack(imagestack{1}, Histology.stack{1}.rec);
    % In case you scanned a second channel uncomment the following line:
    % Histology.stack{1}.Ch2 = crop_stack(imagestackch2{1}, Histology.stack{1}.rec);
else
    % Here, stack 1 and stack 2 are already of the same size. Therefore, we don't need to
    % crop them. Instead, copy imagestack 1 to the final Data structure.
    Histology.stack{1}.Ch{1}=imagestack{1};
end

% Compensate for case of out-of-focus images at the upper or lower borders of the image
% stack. The projection parameters must be estimated. Playback the
% imagestack by implay(Histology.stack{1}.Ch{1})
Stack_begin_MaxIP=max(Histology.stack{1}.Ch{1}(:,:,1:5),[],3);
Stack_end_MaxIP=max(Histology.stack{1}.Ch{1}(:,:,119:end),[],3);
Histology.stack{1}.Ch{1}(:,:,[1:5 119:end])=[];
Histology.stack{1}.Ch{1}=cat(3,Stack_begin_MaxIP, Histology.stack{1}.Ch{1},Stack_end_MaxIP);
clear Stack_begin_MaxIP Stack_end_MaxIP

% reduce image noise
for j=1:size(Histology.stack{1}.Ch{1},3)
    Histology.stack{1}.Ch{1}(:,:,j) = wiener2(Histology.stack{1}.Ch{1}(:,:,j),[3 3]);
end

for j=1:size(Histology.stack{1}.Ch{1},3)
    Histology.stack{1}.Ch{1}(:,:,j) = medfilt2(Histology.stack{1}.Ch{1}(:,:,j),[4 4]);
end

% Calculate projection images of stack 1
Histology.stack{1}.minip.end=min(Histology.stack{1}.Ch{1}(:,:,end-2:end),[],3);
Histology.stack{1}.maxip.end=max(Histology.stack{1}.Ch{1}(:,:,end-2:end),[],3);

% Stack 2
% check if stack 1 and 2 are of different size. If yes - use cropping
% rectangle defined for stack 1
if size(imagestack{1},1)~= size(imagestack{2},1) | size(imagestack{1},2) ~= size(imagestack{2},2)
    % position the rectangle for croppping stack 2
    [t, Histology.stack{2}.rec]=b_croprect(imagestack{2}(:,:,20),Histology.stack{1}.rec);
    % crop stack 2
    stack_crop{2}.Ch{1} = crop_stack(imagestack{2}, Histology.stack{2}.rec);
else
    % Here, stack 1 and stack 2 are already of the same size. Therefore, we don't need to
    % crop them. Instead, copy imagestack 2 to the final Data structure.
    Histology.stack{2}.Ch{1} = imagestack{2};
end

% Compensate for case of out-of-focus images at the upper or lower borders of the image
% stack. The projection parameters must be estimated. Playback the
% imagestack by implay(Histology.stack{2}.Ch{1})
Stack_begin_MaxIP=max(Histology.stack{2}.Ch{1}(:,:,1:20),[],3);
Stack_end_MaxIP=max(Histology.stack{2}.Ch{1}(:,:,end-5:end),[],3);
Histology.stack{2}.Ch{1}(:,:,[1:20 end-5:end])=[];
Histology.stack{2}.Ch{1}=cat(3,Stack_begin_MaxIP, Histology.stack{2}.Ch{1},Stack_end_MaxIP);
clear Stack_begin_MaxIP Stack_end_MaxIP

% reduce image noise
for j=1:size(Histology.stack{2}.Ch{1},3)
    Histology.stack{2}.Ch{1}(:,:,j) = wiener2(Histology.stack{2}.Ch{1}(:,:,j),[3 3]);
end

for j=1:size(Histology.stack{2}.Ch{1},3)
    Histology.stack{2}.Ch{1}(:,:,j) = medfilt2(Histology.stack{2}.Ch{1}(:,:,j),[4 4]);
end

% Calculate projection images of stack 2
Histology.stack{2}.minip.start=min(Histology.stack{2}.Ch{1}(:,:,1:2),[],3);
Histology.stack{2}.maxip.start=max(Histology.stack{2}.Ch{1}(:,:,1:2),[],3);

%% Reconstruct 3D histology model
% Co-register the z-stacks of the two histological sections, scanned with a laser confocal microscope

%initialize image registration method
Options=struct('Similarity','cc','Registration','Both','Penalty',0.0005,'MaxRef',3,'Grid',[],'Spacing',[15 15],'MaskMoving',[],'MaskStatic',[],'Verbose',0,'Points1',[],'Points2',[],'PStrength',[],'Interpolation','Linear','Scaling',[1 1]);

%Co-register section 1/2
HREG{1}.fixed=imadjust(Histology.stack{1}.minip.end);
HREG{1}.moving=imadjust(Histology.stack{2}.minip.start);
HREG{1}.Options = Options;

% Select control points for the registration process. Use pre-defined
% control points if available.
if exist('Sampledata','var')
    HREG{1}.Options.Points1= fliplr(Sampledata.HREG{1}.Options.Points1);
    HREG{1}.Options.Points2= fliplr(Sampledata.HREG{1}.Options.Points2);
    [HREG{1}.Options.Points1,HREG{1}.Options.Points2] = cpselect(HREG{1}.moving,HREG{1}.fixed, HREG{1}.Options.Points1,HREG{1}.Options.Points2,'Wait',true);
else
    cpselect(HREG{1}.moving,HREG{1}.fixed,'Wait',false);
end
HREG{1}.Options.Points1=fliplr(HREG{1}.Options.Points1);
HREG{1}.Options.Points2=fliplr(HREG{1}.Options.Points2);
HREG{1}.Options.PStrength=ones(size(HREG{1}.Options.Points1,1),1);
HREG{1}.Options.PStrength=HREG{1}.Options.PStrength.*0.95;

% Perform registration
[HREG{1}.Ireg,HREG{1}.O_trans,HREG{1}.Spacing,~,~,~] = image_registration(HREG{1}.moving,HREG{1}.fixed ,HREG{1}.Options);
HREG{1}.fixed = im2uint8(HREG{1}.fixed);

% Validate registration result
implay(cat(3,HREG{1}.fixed,HREG{1}.Ireg))

%% Transform 2nd image strack and merge them
for i=1:size(Histology.stack{2}.Ch{1},3)
    Histology.stack{2}.aligned.ch{1}(:,:,i)=im2uint8(bspline_transform(HREG{1}.O_trans,Histology.stack{2}.Ch{1}(:,:,i),HREG{1}.Spacing,3));
end
Histology.stack{2}.aligned.ch{1}=mat2gray(Histology.stack{2}.aligned.ch{1});
Histology.stack{2}.aligned.ch{1} = im2uint8(Histology.stack{2}.aligned.ch{1});
Histology.Ch{1}.stack = cat(3,Histology.stack{1}.Ch{1} , Histology.stack{2}.aligned.ch{1});

% remove individual stacks as we now have the final histology model
Histology = rmfield(Histology,'stack');
clear imagestack

% Playback co-registered histology model
implay(Histology.Ch{1}.stack)
%% Crop 3D histology model to the lesion of the endoscopic lens
if exist('Sampledata','var')
    [~, Histology.rec ] = b_croprect(Histology.Ch{1}.stack(:,:,1),   Sampledata.Histology.rec);
    pause (2)
    close
    Histology.Ch{1}.stack = crop_stack(Histology.Ch{1}.stack,  Histology.rec);
else
    [~, Histology.rec ] = b_croprect(Histology.Ch{1}.stack(:,:,1));
    Histology.Ch{1}.stack = crop_stack(Histology.Ch{1}.stack,  Histology.rec);
end

%% Automatically detect blood vessels in the 3D histology model
[Histology.vessels, Histology.Ch{1}.stackEQ] = detectvessels3D(Histology.Ch{1}.stack, 18);

%% Manually segment cell soma in the 3D histology model
if exist('Sampledata','var')
    Histology.segmentedcells = annotatesegmentedcells(Histology.Ch{1}.stack,Histology.Ch{1}.stackEQ, Sampledata.Histology.segmentedcells);
else
    Histology.segmentedcells = annotatesegmentedcells(Histology.Ch{1}.stack,Histology.Ch{1}.stackEQ);
end
%% align Ca2+ and histo images

% 1) find working distance
if exist('Sampledata','var')
    [Multireg, Histology] =findworkingdistance(Histology,RECDSA,Sampledata.Multireg);
else
    [Multireg, Histology] =findworkingdistance(Histology,RECDSA);
end

% 2) rotation optimization
Multireg = rotation_optimization(Histology.vessels, Multireg);

[Histology, Multireg] = createrotatedstack(Histology,Multireg)

% 3) final non-rigid registration 
Multireg = multimodal_alignment(Histology, Multireg,RECDSA);

%% Align in vivo recorded cells with histology
CaAligned = multimodal_ca_alignment(Allcells, RECDSA, Histology, Multireg)

%% match cells
[Cellmatch, Multireg] = match_multimodal_cells(CaAligned,Histology, Multireg)

%% show identified cells
plot_interactive_identifiedcells(Histology,Cellmatch,Multireg)