function [Histocell,Calcium_cells, cellmatch]= multimodal_cell_match_pairs(HistoImageStack, HistoCells, Calcium_cells, cellsize)
% This function calculates matching scores for cells detected in in vivo
% ca2+ imaging and post hoc histology. 
%
% Inputs:
%          HistoImageStack : a 3D matrix representing images from the 3D
%          histology model.
%
%          HistoCells : a 3D logical matrix representing segmented cells in
%          the histology model.
%
%          Calcium_cells : struct containing in vivo recorded neurons
%
%          cellsize : estimated diameter of the cell soma of in vivo 
%                     recorded neurons
%
%
% Outputs:  
%          Histocell : struct containing parameters of cells identified in 
%                      histology for multimodal cell match         
%                      
%          Calcium_cells : struct containing parameters of cells identified in 
%                      in vivo Ca2+ imaging for multimodal cell match         
%
%          cellmatch : metrics for matching in vivo recorded cells post hoc
%          in histology
%
% Function is written by Philip Anner (2020)

if ~isa(HistoCells,'logical')
   HistoCells = imbinarize(HistoCells); 
end


allhistocells = zeros(size(HistoCells));

C = bwconncomp(HistoCells,26);
full3dcell = cell(1,C.NumObjects);
for i=1:(C.NumObjects)
    clear tmpdendrites
    
    %% segmented cell soma
    tmpstack_soma = zeros(size(HistoCells));
    tmpstack_soma(C.PixelIdxList{i})=1;
     
    [~, ~, Zdim] = ind2sub(size(HistoCells),C.PixelIdxList{i});
    Histocell.start(i) = min(Zdim);
    Histocell.end(i)= max(Zdim);
    Histocell.soma(:,:,i)=max(tmpstack_soma,[],3);
    
    Histocell.image(:,:,i) = max(HistoImageStack(:,:,(Histocell.start(i): Histocell.end(i))),[],3);
    c = regionprops3(tmpstack_soma, 'centroid');
    Histocell.centroid(:,:,i)=round(c.Centroid);

    %% area including dendrites
    tmporigsubstack = HistoImageStack(:,:,Histocell.start(i):Histocell.end(i));
    tmpsomasubstack = tmpstack_soma(:,:,Histocell.start(i):Histocell.end(i));
    tmpdendrites = zeros([size(tmpstack_soma,1) size(tmpstack_soma,2) size(tmporigsubstack,3)]);
    for x=1:size(tmporigsubstack,3)
        tempimg=im2uint8(tmporigsubstack(:,:,x));
        tempimg_filtered= imgaussfilt(tempimg, 2);
        
        deg = [0 45 90];
        tempimg_filtered_detect_h = zeros([size(tempimg_filtered) length(deg)]);
        tempimg_filtered_detect_h_conv = zeros([size(tempimg_filtered) length(deg)]);
        for y=1:length(deg)
            se = strel('line',5, deg(y));
            tempimg_filtered_detect_h(:,:,y) = imtophat(tempimg_filtered,se);
            tempimg_filtered_detect_h_conv(:,:,y)=conv2(tempimg_filtered_detect_h(:,:,y), double(se.getnhood),'same');
        end
        
        
        tI=max(tempimg_filtered_detect_h_conv,[],3);
        tI=mat2gray(tI);
        tI=im2double(tI);
        tI= imgaussfilt(tI, 1);
        tI=imadjust(tI);
        %tI= im2bw(tI, (graythresh(tI)*1.8));
        tI= imbinarize(tI);

        tI=bwareaopen(tI, 100);
        tI=tI-bwareaopen(tI, 1000);
        tmpdendrites(:,:,x)=tI;
    end
    
    fullcell = zeros(size(tmpdendrites));
    Ctmpdendrites = bwconncomp(tmpdendrites ,6);
    Ctmpcellsoma = bwconncomp(tmpstack_soma(:,:,Histocell.start(i):Histocell.end(i)));
    
    clear tmpstack
    founddendrites=0;
    for y=1:Ctmpdendrites.NumObjects
        if size(Ctmpcellsoma.PixelIdxList,2) >1
            error('more than one cell soma detected!')
        end
        
        if  max(ismember(Ctmpcellsoma.PixelIdxList{1},Ctmpdendrites.PixelIdxList{y}) )
            tmpstack = zeros(size(tmpdendrites));
            tmpstack(Ctmpdendrites.PixelIdxList{y})=1;
            fullcell = max(fullcell,tmpstack);
            founddendrites = 1;
        end
        
        if founddendrites
         fullcell= max(tmpsomasubstack, tmpstack);  
        else
            fullcell = tmpsomasubstack;
        end
    end
    Histocell.fullcell(:,:,i) = max(fullcell,[],3);
    t = regionprops3(Histocell.fullcell(:,:,1), 'Orientation');
    Histocell.orientation = t.Orientation;
   
    %% matrix for all cells with dendrites

    pad1 = zeros(size(HistoCells,1),size(HistoCells,2),max(1:Histocell.start(i))-1);
    pad2 = zeros(size(HistoCells,1),size(HistoCells,2),max(size(HistoCells,3)-Histocell.end(i)));
    fullcell = cat(3,pad1,fullcell,pad2);
    fullcell=(fullcell>0);
    fullcell=fullcell*i;
    allhistocells = max(allhistocells, fullcell);
    
    %reshape to 2d matrix and make sparse for memory efficiency
    t=reshape(fullcell,size(allhistocells,1),[],1);
    full3dcell{i}=sparse(t);

    %% PDF calculation - option of region based PDF kernel
    % get pdf per cell
    tI = zeros(size(Histocell.image(:,:,1)));
    tI(Histocell.centroid(:,2,i),Histocell.centroid(:,1,i))=1;
    
    % use gauss kernel size depending on local matching confidence by local
    % jaccard values
    Histocell.centroidimage(:,:,i) = tI;
    Histocell.pdf(:,:,i) = mat2gray(imgaussfilt(tI, cellsize));
end

% remove  merged cells due to region growing   
for i=1:(C.NumObjects)
    for j=1:(C.NumObjects)
        if i~=j
            full3dcell{i} = full3dcell{i} - full3dcell{j};  
            full3dcell{i} = full3dcell{i} > 0;
        end
    end
end

allhistocells=zeros(size(allhistocells));
for i=1:(C.NumObjects)
    full3dcell{i} = full(full3dcell{i});
    full3dcell{i} = reshape(full3dcell{i}, size(allhistocells,1),size(allhistocells,2),size(allhistocells,3));
    allhistocells = max(allhistocells, full3dcell{i});
end

allhistocells=bwareaopen(allhistocells,100);
Histocell.separatedcells = allhistocells;
clear tempimg_filtered_detect tmpstack_soma pad1 pad2 tempimg_filtered_detect_h tempimg_filtered_detect_h_conv fullcell

%% ca2+ cells
if ~isa(Calcium_cells.areas,'logical')
    Calcium_cells.areas=imbinarize(Calcium_cells.areas);
end
if ~isa(Calcium_cells.soma,'logical')
    Calcium_cells.soma=imbinarize(Calcium_cells.soma);
end
% calculate cell properties
for i=1:size(Calcium_cells.areas,3)
    t = regionprops3(Calcium_cells.areas(:,:,1), 'Orientation');
    Calcium_cells.orientation = t.Orientation;
end

%% Define cell Area or Soma for matching
for i=1:size(Calcium_cells.areas,3)
  
    t=regionprops(Calcium_cells.areas(:,:,i), 'centroid');
    if size(t,1) > 1
        s = regionprops(Calcium_cells.areas(:,:,i), 'Area');
        s = struct2cell(s);
        s = [s{:}];
        [~, ind ]= max(s);
        Calcium_cells.cellareacentroid(i,:) = [t(ind).Centroid i];
    else
        Calcium_cells.cellareacentroid(i,:)=[t.Centroid i];
    end
    t=regionprops(Calcium_cells.soma(:,:,i), 'centroid');
    if size(t,1) > 1
        s = regionprops(Calcium_cells.areas(:,:,i), 'Area');
        s = struct2cell(s);
        s = [s{:}];
        [~, ind ]= max(s);
        Calcium_cells.somacentroid(:,:,i) = [t(ind).Centroid i];
    else
        Calcium_cells.somacentroid(:,:,i)=[t.Centroid i];
    end
end

ca2centers = zeros(size(Calcium_cells.areas));
ind=ceil(Calcium_cells.somacentroid);
for i=1:size(Calcium_cells.somacentroid,3)
    ca2centers(ind(1,2,i),ind(1,1,i),ind(1,3,i))=1;
end
ca2centerspdf = zeros(size(ca2centers));

for i=1:size(ca2centers,3)   
    ca2centerspdf(:,:,i)=imgaussfilt(ca2centers(:,:,i),cellsize);
end
ca2centerspdf=mat2gray(ca2centerspdf);

%% PDF for normalization

refsize=size(Histocell.image(:,:,1));
t=zeros((refsize));
t(round(refsize(1)/2), round(refsize(2)/2))=1;
t=imgaussfilt(t,cellsize); % 10
t=mat2gray(t);
maxPDF=trapz(trapz(t.*t));

%%
if ~isa(Histocell.fullcell,'logical')
    Histocell.fullcell = imbinarize(Histocell.fullcell);
end

if ~isa(Histocell.soma,'logical')
    Histocell.soma = imbinarize(Histocell.soma);
end

if ~isa(Calcium_cells.soma,'logical')
    Calcium_cells.soma = imbinarize(Calcium_cells.soma);
end

if ~isa(Calcium_cells.areas,'logical')
    Calcium_cells.areas = imbinarize(Calcium_cells.areas);
end

clear cellmatch
% for all histo cells
for i=1:(C.NumObjects)
    % for all cells from ca2+
    for j=1:size(Calcium_cells.areas,3)

        t =Histocell.pdf(:,:,i).*ca2centerspdf(:,:,j);
        
        % normalize by PDF from specific histo cell area
        cellmatch.pdf(j,i)=trapz(trapz(t))/maxPDF;
        cellmatch.allareajaccard(j,i)=jaccard(Histocell.fullcell(:,:,i),Calcium_cells.areas(:,:,j));
        cellmatch.somajaccard(j,i)=jaccard(Histocell.soma(:,:,i),Calcium_cells.soma(:,:,j));
        cellmatch.euclidean(j,i) = sqrt((Histocell.centroid(1,1,i)-Calcium_cells.somacentroid(1,1,j))^2 + (Histocell.centroid(1,2,i)-Calcium_cells.somacentroid(1,2,j))^2);
        
        cellmatch.zdim(j,i)=Histocell.centroid(:,3,i);        
        A=round(Calcium_cells.somacentroid(1,1,j))-Histocell.centroid(:,1,i);
        B=round(Calcium_cells.somacentroid(1,2,j))-Histocell.centroid(:,2,i);
        
        [theta,rho] = cart2pol(A,B);
        cellmatch.theta(j,i) =theta;
        cellmatch.rho(j,i) =rho;
        cellmatch.U(j,i)=A;
        cellmatch.V(j,i)=B;        
    end
end


end
