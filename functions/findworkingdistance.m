function [Multireg, Histology]= findworkingdistance(Histology,RECDSA, varargin)
% This function lets the user interactively co-register in vivo Ca2+
% imaging and post hoc histology using linear transformations and to 
% determine the in vivo working distance of the miniscope post hoc in 
% histology.
%
% Inputs:
%           Histology : Data structure containing histology imaging data.
%                       Blood vessels detected in histology must be 
%                       included as Histology.vessels. Call 
%                       detectvessels3D(Histostack, maxvesselsize) first.
%
%           RECDSA : Data structure containing in vivo detected blood 
%                    vessels. Acquired by calling InVivo_DetectBloodVessels
%
%           Multireg (optional) : Multimodal registration data acquired by
%                                 findworkingdistance(Histology,RECDSA).
%                                 Stored control points for registration
%                                 can be loaded by calling: 
%                           findworkingdistance(Histology,RECDSA, Multireg)
%
% Outputs: 
%           Multireg : Data structure for pre-processing and registration
%                      control pints for multimodal image alignment between
%                      histology and in vivo acquired Ca2+ imaging data.
%
%           Histology : Data structure containing histology imaging data.
%
% Functionality : 
%                   Interactively identify the working distance of in
%                   vivo Ca2+ imaging post hoc using blood vessels. 
%       
%                     arrow down: next image
%                     arrow up: previous image
%                     arrow right: jump + 10 images
%                     arrow left: jump - 10 images
%                     a: overlay Ca2+ detected blood vessels
%                     s : show merged Ca2+ detected blood vessels
%
% Function is written by Philip Anner (2020)

closedvideoplayer = [];
closedvideoplayer = 0;
if numel(varargin) == 1
    Multireg = varargin{1};
elseif numel(varargin) > 1
    error('Wring number of input arguments! Call findworkingdistance(Histology,RECDSA) or findworkingdistance(Histology,RECDSA, Multireg)')
else
    Multireg = [];
end

Multireg.preprocess.rot90 = 0;
Multireg.preprocess.cropca = 0;
Multireg.preprocess.fliplr=0;
Multireg.preprocess.matchimagesize=0;


if isfield(Histology,'fit') && isfield(Histology.fit,'vessels')
    Hvessels = Histology.fit.vessels;
elseif isfield(Histology,'vessels')
    Hvessels = Histology.vessels;
else
    error('Blood vessels not found in Histology. Please detect blood vessels by calling detectvessels3D(Histostack, maxvesselsize)!')
end

if isfield(RECDSA,'vesselimg_aln')
    cavessels = RECDSA.vesselimg_aln;
elseif isfield(RECDSA,'vesselimg')
    cavessels = RECDSA.vesselimg;
else
    error('Blood vessels detected in Ca2+ imaging not found. Please provide blood vessels in RECDSA.vesselimg by calling InVivo_DetectBloodVessels')
end

answer = questdlg('Would you like to do an initial alignment?', ...
    'Yes', ...
    'No');
% Handle response
switch answer
    case 'Yes'
        camod = 1;
        histomod = 1;
        F = figure();
        ax1 = subplot(1,2,1);
        imshow(cavessels);
        title({'Ca2+ - vessels', ['dimensions: ' num2str(size(cavessels,1)) ' x ' num2str(size(cavessels,2))]})
        ax2 = subplot(1,2,2);
        imshow(max(Hvessels(:,:,1:round(size(Hvessels,3)/3)),[],3));
        title({'2D projection of blood vessels detected', ['in Histology using image 1 to ' num2str(round(size(Hvessels,3)/3))], ['dimensions: ' num2str(size(Hvessels,1)) ' x ' num2str(size(Hvessels,2))]})
        msgbox('Press any key to continue')
        pause
        %waitforbuttonpress
        
        while camod
            answerca = questdlg('Ca2 image - transformations?','Configuration','rot90', 'fliplr', 'done','done');
            switch answerca
                case 'rot90'
                    Multireg.preprocess.rot90 = Multireg.preprocess.rot90 +1;
                    cavessels = rot90(cavessels);
                    axes(ax1), imshow(cavessels)
                case 'fliplr'
                    if Multireg.preprocess.fliplr==0
                        Multireg.preprocess.fliplr = 1;
                    else
                        Multireg.preprocess.fliplr = 0;
                    end
                    
                    cavessels = fliplr(cavessels);
                    axes(ax1), imshow(cavessels)
                    % case 'crop'
                    %  cavessels = b_croprect(cavessels, [1 size(Hvessels,2) 1 size(Hvessels,1)]);
                case'done'
                    camod = 0;
            end
        end
        
        while histomod
            answerhisto = questdlg('Histo image - transformations?','Configuration','crop','done','done');
            switch answerhisto
                case 'crop'
                    if size(Hvessels(:,:,1)) ~= size(cavessels)
                        [~, recH]= b_croprect(max(Hvessels(:,:,1:round(size(Hvessels,3)/3)),[],3), [1 1 size(cavessels,2) size(cavessels,1)]);
                        Hvessels= imcrop3D(Hvessels, recH);
                        
                        ax2s(ax1), imshow(Hvessels(:,:,1));
                        axes(ax2), imshow(cavessels);

                        implay(Hvessels)
                        Multireg.preprocess.crophisto = 1;
                        Histology.fit.vessels = Hvessels;
                       
                        if isfield (Histology,'Ch')
                            for i = 1:size(Histology.Ch,2)
                                Histology.fit.Ch{i} =  imcrop3D(Histology.Ch{i}.stack,recH); % crop stack
                                Histology.fit.Ch{i} =  imcrop3D(Histology.Ch{i}.stackEQ,recH); % crop stack
                            end
                        end
                        
                        if isfield (Histology,'ch')
                            for i = 1:size(Histology.ch,2)
                                Histology.fit.ch{i} =  imcrop3D(Histology.ch{i}.stack,recH); % crop stack
                                Histology.fit.ch{i} =  imcrop3D(Histology.ch{i}.stackEQ,recH); % crop stack
                            end
                        end
                        Histology.histocroprec = recH;
                    else
                        msgbox('images have already same size')
                    end
                    
                case 'done'
                    histomod = 0;
            end
        end
        
        close(F)
        Hvessels_t = (max(Hvessels(:,:,1:round(size(Hvessels,3)/3)),[],3));
        
        
        if size(Hvessels_t) ~= size(RECDSA.vesselimg)
            [ ~ , cavessels ] = matchimagesize(Histology.vessels(:,:,1),cavessels);
            Multireg.preprocess.matchimagesize = 1;
        end
        
        
        regmod = 1;
        while regmod
            
            if isfield(Multireg.preprocess,'fixedPoints')
                [movingPoints, fixedPoints] = cpselect(cavessels,Hvessels_t,Multireg.preprocess.movingPoints ,Multireg.preprocess.fixedPoints, 'Wait',true);
                
            elseif exist('movingPoints','var')
                [movingPoints, fixedPoints] = cpselect(cavessels,Hvessels_t, movingPoints,fixedPoints, 'Wait',true);
            else
                [movingPoints, fixedPoints] = cpselect(cavessels,Hvessels_t,'Wait',true);
            end
            
            Multireg.preprocess.movingPoints = movingPoints;
            Multireg.preprocess.fixedPoints = fixedPoints;
            
            tform = fitgeotrans(movingPoints,fixedPoints,'affine');
            ref = imref2d(size(Hvessels_t));
            cavessels_registered = imwarp(cavessels,tform,'OutputView',ref);
            
            
            h = implay(cat(3,Hvessels_t,cavessels_registered));
            F = h.Parent;
            msgbox('Validate registration result. Close figure to continue');
            waitfor(F);
         
            answerreg = questdlg('Would you like to do the alignment again?', ...
                'Yes', ...
                'No');
            
            switch answerreg
                case 'No'
                    regmod = 0;
            end
        end
    case 'No'
        cavessels_registered = cavessels;
end





steps = 1:size(Histology.vessels,3);
h = waitbar(0,'Calculating projection images...');
stepimages = zeros([size(Hvessels(:,:,1)) size(steps,2)]);
for i=1:size(steps,2)
    stepimages(:,:,i) = (max(Hvessels(:,:,1:steps(i)),[],3));
    waitbar(i/size(steps,2));
end
close(h)

e = 1;
switchvessels = 0;
overlaymode = 0;
i=1;
F = figure;
m = msgbox({'Please specify working distance.', 'Controls: ', 'arrow down: next image',  'arrow up: previous image','arrow right: jump + 10 images', 'arrow left: jump - 10 images', 'a: overlay Ca2+ detected blood vessels', 's : show merged Ca2+ detected blood vessels'});
waitfor(m)
figure(F);
while e
    
    if ~switchvessels
        imshow(stepimages(:,:,i));
    else
        imshow(cavessels_registered);
    end
    
    if overlaymode
        imshowpair(stepimages(:,:,i),cavessels_registered);
    end
    
    %viscircles to prevent figure resize
    viscircles(([0 0 ]), 1,'EdgeColor','b');
    title(['MAX(1:' num2str(steps(i)) ') - step Nr: ' num2str(i) '. To accept selected depth - press Q'])
    
    waitforbuttonpress;
    h=gcf;
    val=double(get(h,'CurrentCharacter'));
    
    try
        switch val
            %upArrow
            case 30
                i=i-1;
                if i == 0
                    i=i+1;
                end
                
                %downArrow
            case 31
                i=i+1;
                if i > size(steps,2)
                    i=i - 1;
                end
                
            %right arrow
            case 29
                i=i+10;
                if i > size(steps,2)
                    i=i - 10;
                end
            %left arrow    
            case 28
                 i=i-10;
                if i <= 0
                    i=i+1;
                end
                
                
                % a - show calcium vessels
            case 97
                if ~switchvessels
                    switchvessels = 1;
                else
                    switchvessels = 0;
                end
                % s - overlay mode
            case 115
                if ~overlaymode
                    overlaymode = 1;
                else
                    overlaymode = 0;
                end
                
            case 113
                e = 0;
                otherwisey
                sprintf('invalid command')
        end
    catch
        warning('unknown command');
    end
end

Multireg.cavesselsAffinereg = cavessels_registered;
Multireg.Histo.vesselimgdepth = (steps(i));
close

 function closecallback(src,callbackdata)
                closedvideoplayer = 1;
 end

end