% Process final tail vein injection experiment
%
% Load projection images, exported from Inscopix Mosaic software
%
% Inat_std : native (prior to fluorescence dye injection) standard deviation 
% projection image.
%
% ICE_std : contrast enhanced (following fluorescence dye injection) 
% standard deviation projection image.
%
% Inat_min : native (prior to fluorescence dye injection) minimum
% intensity projection image.
%
Inat_std = imread('.\Data\In vivo Ca2+\1912-DSA\JP68-1912-mean-native.tif'); 
ICE_std = imread('.\Data\In vivo Ca2+\1912-DSA\JP68-1912-MaxIP_CE.tif'); 
Inat_min = imread('.\Data\In vivo Ca2+\1912-DSA\JP68-1912-min_native.tif'); 

Inat_min = mat2gray(Inat_min);
Inat_std=imgaussfilt(Inat_std,4);
ICE_std=imgaussfilt(ICE_std,4);

DSA=Inat_std-ICE_std;
figure, imshow(DSA,[])
title('Difference image (in vivo)')
DSA=mat2gray(DSA);
pause(2)
close

% topological filtering
clear tempimg_filtered_detect
vesselsenhanced = zeros(size(DSA));
% You may want to increase strength if you have large blood vessels in the
% field of view.
strength = (1:2:15);

% Topological filtering operation to enhance tubular blood vessel
% structures
    for k=1:length(strength)
        deg = [0 25 45 75 90 115 135 150 180 ];
        for j=1:length(deg)
            se = strel('line',strength(k), deg(j));
            tempimg = imbothat(DSA,se);
            vesselsenhanced = max(tempimg,vesselsenhanced);
        end
    end
    
vesselsenhanced= imadjust(vesselsenhanced);

% Use Hessian to detect Blood vessel structures
options = struct('FrangiScaleRange', [1 8], 'FrangiScaleRatio', 2, 'FrangiAlpha', 0.0000000000001, 'FrangiBeta', 1000, 'FrangiC', 50, 'verbose',true,'BlackWhite',false);
vesselsCA=FrangiFilter2D(vesselsenhanced,options);
vesselsCA=mat2gray(vesselsCA);
vesselsCA= imadjust(vesselsCA,[0 0.8], [0 1 ],0.52 );

% show detected blood vessels
figure,imshow(vesselsCA,[])
title('Enhanced blood vessels (in vivo)')
pause(2)
close

% binary segment blood vessels
vesselsCAbin=imbinarize(vesselsCA);
vesselsCAbin= bwareaopen(vesselsCAbin, 50);
figure,imshow(vesselsCAbin,[])
title('Segmented blood vessels (in vivo)')
pause(2)
close
%% include in merged sessions
RECDSA.minrec  = Inat_min;
RECDSA.vesselimg = vesselsCA;
RECDSA.DSA  = DSA;
RECDSA.vesselimg_bw =vesselsCAbin;


%% align images
if exist('Sampledata','var')
    % For the specific sample data images from the tail vein injection must be
    % rotated, flipped, and cropped.
    RECDSA.minrec=rot90(RECDSA.minrec);
    RECDSA.minrec=fliplr(RECDSA.minrec);
    
    RECDSA.vesselimg=rot90(RECDSA.vesselimg);
    RECDSA.vesselimg=fliplr(RECDSA.vesselimg);
    
    RECDSA.DSA=rot90(RECDSA.DSA);
    RECDSA.DSA=fliplr(RECDSA.DSA);
    
    RECDSA.vesselimg_bw=rot90(RECDSA.vesselimg_bw);
    RECDSA.vesselimg_bw=fliplr(RECDSA.vesselimg_bw);
    
    % crop images
    [~, RECDSA.preprocess.rect]= b_croprect(RECDSA.minrec, Sampledata.RECDSA.rec);
    RECDSA.preprocess.cropca = 1;
    RECDSA = registerrecording(RECDSA,RECDSA);
    pause(2)
    close
end
clear Inat_min vesselsCA DSA vesselsCAbin Inat_std P  strength deg options j k vesselsenhanced se tempimg ICE_std
