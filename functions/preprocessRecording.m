function [REC] = preprocessRecording(folder, unmixing, events)
% This function reads in vivo ca2+ imaging data, processed by Inscopix 
% Mosaic software.
%
% Inputs:
%           folder : The folder containing a .csv file with the calcium
%           transients. The .csv file can be exported from Inscopix Mosaic 
%           and must contain variable names (line 1), a time code (column 1),
%           two two columns for each cell containing a dF/F trace, and
%           a Calcium event trace.
%
%           unmixing : must be 0 or 1 if you want to use un-mixing images
%           exported from Inscopix Mosaic software
%
%           events : specify if the recording contains only dF/F traces 
%                    (events = 0)or also calcium events (events = 1).
%
%
%
%
% Function is written by Philip Anner (2020)

% set up variables
convertICmfromdouble = false;
convertICufromdouble = false;

if ~strcmp(folder(end), '\')
    folder = strcat(folder,filesep);
end

tracepath=folder;

t=dir([tracepath '*traces.csv']);

if ~isempty(t)
    f=fopen([tracepath t.name]);
    tmptext=textscan(f, '%s', 1, 'delimiter', '\n');
    fclose(f);
    tt=cell2mat(tmptext{1});
    cellnames=textscan(tt,'%s', 'delimiter', ',');
else
    error(['*.traces file not found at path: ' tracepath ])
end

if events
    % extract cell names
    IC_names = cell((size(cellnames{1},1)-1)/2,1);
    ind=1;
    for i=2:2:size(cellnames{1},1)
        IC_names(ind) =cellnames{1}(i,1);
        ind=ind+1;
    end
    % ca2 traces
    tmptraces = csvread([tracepath t.name],1);
    traces.timing = tmptraces(:,1);
    ind=1;
    for i=2:2:size(cellnames{1},1)
        traces.trace{ind} = tmptraces(:,i);
        traces.events{ind} = tmptraces(:,i+1);
        ind=ind+1;
    end
else
    IC_names = cell((size(cellnames{1},1)-1),1);
    ind=1;
    for i=2:size(cellnames{1},1)
        IC_names(ind) =cellnames{1}(i,1);
        ind=ind+1;
    end
    % ca2 traces
    tmptraces = csvread([tracepath t.name],1);
    traces.timing = tmptraces(:,1);
    ind=1;
    for i=2:size(cellnames{1},1)
        traces.trace{ind} = tmptraces(:,i);
        ind=ind+1;
    end
end
    


%index for images containing ICs
ind = 1;

if ~strcmp(folder(end), '\')
     folder = strcat(folder,filesep);
end
files = dir([folder '*IC*']);
dirFlags = [files.isdir];
ICfolder = files(dirFlags);

if size(ICfolder,1)>1 || size(ICfolder,1)==0 
   F = uigetdir(folder,'Non or multiple IC folders found! Please select IC folder');
   % extract folder for recording data and IC folder
   [s,e] = regexp(F,'.*\\');
   clear ICfolder
   ICfolder.folder = F(s:e);
   ICfolder.name = F(e+1:end);
end


if ~isempty(ICfolder)
    s = dir([ICfolder.folder '\' ICfolder.name '\*IC*.tif']);
    if ~isempty(s)
        filelist = {s.name}';
        filelist = sort_nat(filelist);
    else
        error(['IC images not found at path: ' tracepath ])
    end
else
    error(['IC folder not found at: ' folder ])
end

%% check for minrec image
clear t
t = dir([tracepath '*minrec.tif']);

if size(t,1) == 0
    error('minrec not found!')
end

for i=1:size(t,1)
    if regexp(t(i).name, '.*dff.*')
        IC_dff_minrec = mat2gray(imread([tracepath t(i).name]));
    else
          minrec=mat2gray(imread([tracepath t(i).name]));
          IC_minrec = minrec;
    end
end
    
%% check for maxrec image
clear t
t = dir([tracepath '*maxrec.tif']);
if size(t,1) > 0
    for i=1:size(t,1)
        if regexp(t(i).name, '.*dff.*')
            IC_dff_maxrec = mat2gray(imread([tracepath t(i).name]));
        else
            maxrec=mat2gray(imread([tracepath t(i).name]));
            IC_maxrec = maxrec;
        end
    end
end
    
%% check for meanrec image
clear t
t = dir([tracepath '*meanrec.tif']);

if size(t,1) > 0
    for i=1:size(t,1)
        if regexp(t(i).name, '.*dff.*')
            IC_dff_meanrec = mat2gray(imread([tracepath t(i).name]));
        else
            meanrec=mat2gray(imread([tracepath t(i).name]));
            IC_meanrec = meanrec;
        end
    end
end
%% check for stdrec image
clear t
t = dir([tracepath '*stdrec.tif']);
if size(t,1) > 0
    for i=1:size(t,1)
        if regexp(t(i).name, '.*dff.*')
            IC_dff_stdrec = mat2gray(imread([tracepath t(i).name]));
        else
            stdrec=mat2gray(imread([tracepath t(i).name]));
            IC_stdrec = stdrec;
        end
    end
end
%%
ii=1;
if unmixing == 1 
    for i=1:2:length(filelist)
        ICm(:,:,ii) = imread( [ICfolder.folder '\' ICfolder.name '\' cell2mat(filelist(i))]);
        ICu(:,:,ii) = imread( [ ICfolder.folder '\' ICfolder.name '\' cell2mat(filelist(i+1))]);
        ii=ii+1; 
    end
else 
     for i=1:length(filelist)
        ICm(:,:,ii) = imread( [ ICfolder.folder '\' ICfolder.name '\' cell2mat(filelist(i))]);
        ICu(:,:,ii) = imread( [ ICfolder.folder '\' ICfolder.name '\' cell2mat(filelist(i))]);
        ii=ii+1;
     end
end

if isa(ICm,'single') 
    convertICmfromdouble = true;
end

if isa(ICu,'single') 
    convertICufromdouble = true;
end

if isa(ICm,'uint8')
    ICm=im2double(mat2gray(ICm));
end
if isa(ICu,'uint8')
    ICu=im2double(mat2gray(ICu));
end

% turn negative values to 0
for i=1:size(ICm,3)
    t = ICm(:,:,i);
    t(t<0)=0;
    if convertICmfromdouble
        t=im2double(mat2gray(t));
    end
    ICm(:,:,i)=t;
end

for i=1:size(ICu,3)
    t = ICu(:,:,i);
    t(t<0)=0;
    if convertICufromdouble
        t=im2double(mat2gray(t));
    end
    ICu(:,:,i)=t;
end

disp([ num2str(ii-1) ' ICs found!!']);

% Allocate memory for variables
ICm_bw = zeros(size(ICm));
ICu_bw = zeros(size(ICm));
ICm_soma = zeros(size(ICm));
ICu_soma = zeros(size(ICm));
IC_rot = zeros(size(ICm));
IC_bw = zeros(size(ICm));
IC_soma = zeros(size(ICm));
IC_perim = zeros(size(ICm));
IC_areas = zeros(size(ICm));
IC_overlay = zeros([size(ICm,1) size(ICm,2) 3 size(ICm,3)]);

for i=1:size(ICm,3)   
   
    % adapt intensity values - mixing image
    tI=imgaussfilt(ICm(:,:,i),4);
    tI=uint8(tI);
    ICm_bw(:,:,i) = imbinarize((tI));
    % convert to BW and morphological image opening to detect soma
    ICm_bw(:,:,i) = imfill(ICm_bw(:,:,i),'holes');
    ICm_bw(:,:,i) = imopen(ICm_bw(:,:,i), ones(5,5));
    ICm_bw(:,:,i) = bwareaopen(ICm_bw(:,:,i), 60);
    ICm_soma(:,:,i)=ICm_bw(:,:,i);
    
     % adapt intensity values - unmixing image
    tI=imgaussfilt(ICu(:,:,i),4);

    ICu_bw(:,:,i) = imbinarize((tI));
    
    % convert to BW and morphological image opening to detect soma
    ICu_bw(:,:,i) = imfill(ICu_bw(:,:,i),'holes');
    ICu_bw(:,:,i) = imopen(ICu_bw(:,:,i), ones(5,5));
    ICu_bw(:,:,i) = bwareaopen(ICu_bw(:,:,i), 60);
    ICu_soma(:,:,i)=ICu_bw(:,:,i);
    
    Am =regionprops(ICm_soma(:,:,i),'Area');
    Au = regionprops(ICu_soma(:,:,i),'Area');

    
    % check if mixed or unmixed image have more than 1 center detected
    if size(Am,1) >1 && size(Au,1)==1
        Am=[];
    elseif size(Au,1) >1 && size(Am,1)==1
        Au=[];
    end
       
    if size(Am,1)>0
        Am = Am.Area;
    else
        Am=0;
    end
        
    if size(Au,1)>0
        Au = Au.Area;
    else
        Au=0;
    end
    
    % compare soma sizes, of IC - use image with smaller soma size,
    % assuming better image quality
    
ok = 1;
    if (Am == 0 && Au== 0)
        IC_rot(:,:,i) = zeros(size(ICm(:,:,i)));
        IC_bw(:,:,i)=  zeros(size(ICm(:,:,i)));
        IC_soma(:,:,i)=  zeros(size(ICm(:,:,i)));
        ok=0;
    end
    
    if  (Am > 1 && Au == 0)
        IC_rot(:,:,i)=ICm(:,:,i);
        IC_bw(:,:,i)=ICm_bw(:,:,i);
        IC_soma(:,:,i)= ICm_soma(:,:,i);
         ok=0;
    end
    
    if  (Am == 0 && Au > 1)
        IC_rot(:,:,i)=ICu(:,:,i);
        IC_bw(:,:,i)=ICu_bw(:,:,i);
        IC_soma(:,:,i)= ICu_soma(:,:,i);
         ok=0;
    end

if ok==1
    if   (Am <= Au && Am > 0)
        IC_rot(:,:,i)=ICm(:,:,i);
        IC_bw(:,:,i)=ICm_bw(:,:,i);
        IC_soma(:,:,i)= ICm_soma(:,:,i);
        disp('MIX')
    elseif Au > 0
        IC_rot(:,:,i)=ICu(:,:,i);
        IC_bw(:,:,i)=ICu_bw(:,:,i);
        IC_soma(:,:,i)= ICu_soma(:,:,i);
        disp('UNMIX')
    end
end

    % detect dendrites
    tempimg=IC_rot(:,:,i);
    tempimg_filtered= imgaussfilt(tempimg, 2);
    
    
    count = 1;
    strength = [ 2 3 4 5 6 7 8 10 15];
    deg = [0 45 90 125 180 225 270 315 300 ];
    tempimg_filtered_detect = zeros([size(tempimg_filtered) numel(strength)* numel(deg)]);
    for k=1:length(strength)
        for j=1:length(deg)
            se = strel('line',strength(k), deg(j));
            tempimg_filtered_detect(:,:,count) = imtophat(tempimg_filtered,se);
            tempimg_filtered_detect(:,:,count)=conv2(double(tempimg_filtered_detect(:,:,count)), double(se.getnhood)','same');
            count = count + 1;
        end
    end
    
    tI=max(tempimg_filtered_detect,[],3);
    tI=mat2gray(tI);
    tI=im2double(tI);
    tI= imgaussfilt(tI, 1);
    tI=imadjust(tI);
    
    mergedI=max(IC_bw(:,:,i),tI);    
    c=regionprops(IC_bw(:,:,i),'Centroid');

    if ~isempty(c)
        c=c.Centroid;
        t=graythresh(double(mergedI));
        [~, J] = regionGrowing(double(mergedI), [round(c(2)) round(c(1))], t, floor(min(size(mergedI))/4), true);
        % if the area of a segmented cell is equal to or greater than the
        % image size, it is likely that region growing segmentes background
        % signals. Modify threshold
        if sum(sum(J>0)) >= numel(mergedI)/2
            [~, J] = regionGrowing(double(mergedI), [round(c(2)) round(c(1))], t/2, Inf, true);
        end
        
        IC_perim(:,:,ind) = bwperim(J);
        IC_areas(:,:,ind) = (J);
        IC_overlay(:,:,:,ind) = imoverlay(IC_rot(:,:,i), IC_perim(:,:,ind), [.3 1 .3]);
        
    else
        IC_perim(:,:,ind) = zeros(size(tI));
        IC_areas(:,:,ind) = zeros(size(tI));
        IC_overlay(:,:,:,ind)  = imoverlay(zeros(size(tI)),zeros(size(tI)));
        
    end
    ind = ind+1;
    
end
IC_overlay = uint8(IC_overlay);
IC_overlay_movie=immovie(IC_overlay);

% make all IC image
REC.allIC = max(IC_rot,[],3);

%
REC.IC = IC_rot;
REC.perim = IC_perim;
REC.areas  = IC_areas;
REC.overlay_movie =  IC_overlay_movie;
REC.soma = IC_soma;
REC.names = IC_names;
REC.traces = traces;
if exist('IC_minrec','var')
    REC.minrec = IC_minrec;
end
if exist('IC_maxrec','var')
    REC.maxrec = IC_maxrec;
end
if exist('IC_meanrec','var')
    REC.meanrec = IC_meanrec;
end
if exist('IC_stdrec','var')
    REC.stdrec = IC_stdrec;
end
if exist('IC_dff_minrec','var')
    REC.dff.minrec = IC_dff_minrec;
end    
if exist('IC_dff_maxrec','var')
    REC.dff.maxrec = IC_dff_maxrec;
end    
if exist('IC_dff_meanrec','var')
    REC.dff.meanrec = IC_dff_meanrec;
end    
if exist('IC_dff_stdrec','var')
    REC.dff.stdrec = IC_dff_stdrec;
end    
end