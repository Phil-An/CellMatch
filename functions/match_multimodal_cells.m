function [output, Multireg] = match_multimodal_cells(Ca2cells, Histology, Multireg)
% This function automatically matches in vivo recorded neurons post hoc in
% Histo imaging data. 
%
% Inputs:
%           Ca2cells :  Anligned in vivo Ca2+ imaging data. Run 
%           multimodal_ca_alignment(Allcells, RECDSA , Histo, Multireg) 
%           for obtaining aligned Ca2cells.
%
%           Histology : Data structure containing Histo imaging data.
%                       
%           Multireg : Data structure containing parameters for multimodal 
%                      image alignment between Histo and in vivo 
%                      acquired Ca2+ imaging data.
%
% Outputs:
%           output : data structure containing variables of matched cell
%           pairs
%           
%          Multireg : Data structure containing multimodal registration 
%                     parameters. 
%        
% Function is written by Philip Anner (2020)

if Multireg.viewpoint.values.rotx ~= 0 || Multireg.viewpoint.values.roty ~= 0 
    if isfield(Histology,'rot')
        Histo = Histology.rot;
    else
        error('Rotation optimization used but rotated histology structure found! Run createrotatedstack(Histology, Multireg) first');
    end
else
    Histo = Histology;
end
clear Histology;

if ~isa(Ca2cells.soma,'logical')
    Ca2cells.soma = imbinarize(Ca2cells.soma);
end

if ~isa(Ca2cells.areas,'logical')
    Ca2cells.areas = imbinarize(Ca2cells.areas);
end

%% Calculate average Ca2+ cell diameter
cellsize = zeros(1,size(Ca2cells.areas,3));
for i=1:size(Ca2cells.areas,3)
    
    if sum(sum(Ca2cells.soma(:,:,i))) > 0
        c = regionprops(Ca2cells.soma(:,:,i), 'Eccentricity');
        Ca2cells.eccentricity(i) = c.Eccentricity;
        
        if Ca2cells.eccentricity(i) > 0.5
            t= regionprops(Ca2cells.soma(:,:,i) , 'MinorAxisLength');
            cellsize(i) = t.MinorAxisLength;
        else
            t= regionprops(Ca2cells.soma(:,:,i) , 'EquivDiameter');
            cellsize(i) = t.EquivDiameter;
        end
    else
        Ca2cells.eccentricity(i)  =NaN;
        cellsize(i) = NaN;
    end
end

Ca2cells.Meancellradius = nanmean(cellsize/2);
sprintf(['Mean cell radius: ' int2str(floor(Ca2cells.Meancellradius)) ' px.'])
clear cellsize

%% calculate distances of cells in ca2+ imaging
for i=1:size(Ca2cells.areas,3)
    for j=1:size(Ca2cells.areas,3)
        Ca2cells.somaeucD(i,j) = sqrt((Ca2cells.somacentroid(1,1,i)-Ca2cells.somacentroid(1,1,j))^2 + (Ca2cells.somacentroid(1,2,i)-Ca2cells.somacentroid(1,2,j))^2);
    end
end
Ca2cells.somaeucD(Ca2cells.somaeucD==0)=NaN;

%% create mask for Histo images restricting to ROI

%% Generate Mask of endoscope ROI
if isfield(Histo,'Ch')
    if isfield(Histo.Ch{1},'stack')
        stack = Histo.Ch{1}.stack;
    end
    
elseif isfield(Histo,'ch')
    if isfield(Histo.ch{1},'stack')
        stack = Histo.ch1.stack;
    end
else
    error('Histology channel 1 not found!')
end

waitfor(msgbox('Please define ROI of endoscope'))

figure, imshowpair(max(stack(:,:,1:30),[],3),max(Ca2cells.soma,[],3))
vesselimgdepth = size(stack,3);
title('Please draw ROI of endoscope. Double click in ROI to finish drawing.')
hold on
h=imellipse;
wait(h); % get position
pos=h.getPosition;
Histo.maskroi=h.createMask;
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','r', 'LineWidth',3);
pause(2);
close
Histo.maskroi=double(Histo.maskroi);
Histo.maskroi(Histo.maskroi==0)=nan;

%% skeletonize blood vessel images
if ~isa(Ca2cells.RECDSA.vesselimg_bw,'logical')
    Ca2cells.RECDSA.vesselimg_bw = imbinarize(Ca2cells.RECDSA.vesselimg_bw,'adaptive');
end
Ca2cells.vesselsCAskel = bwskel(Ca2cells.RECDSA.vesselimg_bw);

Histo.vesselimg = max(Histo.vessels(:,:,1:vesselimgdepth),[],3);
Histo.vesselimg=imadjust(Histo.vesselimg);
Histo.vesselimg=mat2gray(Histo.vesselimg);
Histo.vesselimg_bw = imbinarize(Histo.vesselimg,'adaptive');
Histo.vesselimg_bw = bwareaopen(Histo.vesselimg_bw,200);
Histo.vesselimgbwskel = bwskel(Histo.vesselimg_bw);
Histo.vesselimgbwskel= Histo.vesselimgbwskel .*Histo.maskroi;

% calculate directional filters
[~,~,H.vert,H.hor] = edge(Histo.vesselimgbwskel);
H.vert = imbinarize(H.vert);
H.vert = bwareaopen(H.vert,10);

H.hor = imbinarize(H.hor);
H.hor = bwareaopen(H.hor,10);
% figure, imshow((H.vert))
% figure, imshow(H.hor)


H.v.comp = bwconncomp(H.vert);
H.h.comp = bwconncomp(H.hor);


H.v.orientation = regionprops(H.vert,'Orientation');
H.h.orientation = regionprops(H.hor,'Orientation');

[~,~,Ca.vert,Ca.hor] = edge(Ca2cells.vesselsCAskel);
Ca.vert = imbinarize(Ca.vert);
Ca.vert = bwareaopen(Ca.vert,10);

Ca.hor = imbinarize(Ca.hor);
Ca.hor = bwareaopen(Ca.hor,10);

Ca.v.comp = bwconncomp(Ca.vert);
Ca.h.comp = bwconncomp(Ca.hor);
Ca.v.orientation = regionprops(Ca.vert,'Orientation');
Ca.h.orientation = regionprops(Ca.hor,'Orientation');

%% calculate distance transforms
Ca.vertD = bwdist(Ca.vert,'euclidean');
Ca.horD = bwdist(Ca.hor,'euclidean');
Ca.vertD=imcomplement(Ca.vertD);
Ca.horD=imcomplement(Ca.horD);

Ca.selfD = bwdist(Ca2cells.vesselsCAskel,'euclidean');
Ca.selfD=imcomplement(Ca.selfD);
Ca.comp = bwconncomp(Ca2cells.vesselsCAskel);

H.vertD = bwdist(H.vert,'cityblock');
H.horD = bwdist(H.hor,'cityblock');
H.vertD=imcomplement(H.vertD);
H.horD=imcomplement(H.horD);


[Ca.sGx ,Ca.sGy] = imgradientxy(Ca.selfD);

[Ca.vertGx ,Ca.vertGy] = imgradientxy(Ca.vertD);
[Ca.horGx, Ca.horGy] = imgradientxy(Ca.horD);

[H.vertGx ,H.vertGy] = imgradientxy(H.vertD);
[H.horGx ,H.horGy] = imgradientxy(H.horD);

clear t
for cellind = 1:size(Ca2cells.soma,3)
    
    
    cellcenter = [round(Ca2cells.somacentroid(:,1,cellind)),round(Ca2cells.somacentroid(:,2,cellind))];
    
    
    % find closest vessel in calcium imaging - vertical
    if cellcenter(2) >  size(Ca.vertGx,1)
        cellcenter(2) =  size(Ca.vertGx,1);
    end
    
    if cellcenter(1)> size(Ca.vertGx,2)
        cellcenter(1)= size(Ca.vertGx,2);
    end
    
    Ca.v.posx=[];
    Ca.v.posy=[];
    ind = 1;
    Ca.v.posx(ind)= cellcenter(1);
    Ca.v.posy(ind)= cellcenter(2);
    e=1;
    while e
        if Ca.vertD(Ca.v.posy(ind), Ca.v.posx(ind)) == 1 || Ca.vertD(Ca.v.posy(ind), Ca.v.posx(ind)) == 0
            e=0;
            break
        end
        
        if  Ca.vertGx(Ca.v.posy(ind), Ca.v.posx(ind)) == 0
            Ca.v.posx(ind+1) = Ca.v.posx(ind);
        else
            Ca.v.posx(ind+1) = Ca.v.posx(ind) + round( Ca.vertGx(Ca.v.posy(ind), Ca.v.posx(ind)) /norm(Ca.vertGx(Ca.v.posy(ind), Ca.v.posx(ind))));
        end
        if Ca.vertGy(Ca.v.posy(ind), Ca.v.posx(ind)) == 0
            Ca.v.posy(ind+1) = Ca.v.posy(ind);
        else
            Ca.v.posy(ind+1) = Ca.v.posy(ind) + round( Ca.vertGy(Ca.v.posy(ind), Ca.v.posx(ind)) / norm(Ca.vertGy(Ca.v.posy(ind), Ca.v.posx(ind))));
        end
        
        
        if ind > 1 && (Ca.v.posx(ind)==Ca.v.posx(ind+1) & Ca.v.posy(ind) ==Ca.v.posy(ind+1)) |  (Ca.v.posx(ind-1)==Ca.v.posx(ind+1) & Ca.v.posy(ind) ==Ca.v.posy(ind+1)) | (Ca.v.posx(ind)==Ca.v.posx(ind+1) & Ca.v.posy(ind-1) ==Ca.v.posy(ind+1))
            Ca.v.posx(ind+1) = Ca.v.posx(ind)+randi([-1,1],1);
            Ca.v.posy(ind+1) = Ca.v.posy(ind)+randi([-1,1],1);
        end        
        ind = ind +1;     
    end
    Ca.v.closestvesselx = Ca.v.posx(end);
    Ca.v.closestvessely = Ca.v.posy(end);
       
    % find closest vessel in calcium imaging - horizontal
    Ca.h.posx=[];
    Ca.h.posy=[];
    ind = 1;
    Ca.h.posx(ind)= cellcenter(1);
    Ca.h.posy(ind)= cellcenter(2);
    e=1;
    while e
        if Ca.horD(Ca.h.posy(ind), Ca.h.posx(ind)) == 1  || Ca.horD(Ca.h.posy(ind), Ca.h.posx(ind)) == 0
            e=0;
            break
        end
        
        if  Ca.horGx(Ca.h.posy(ind), Ca.h.posx(ind)) == 0
            Ca.h.posx(ind+1) = Ca.h.posx(ind);
        else
            Ca.h.posx(ind+1) = Ca.h.posx(ind) + round( Ca.horGx(Ca.h.posy(ind), Ca.h.posx(ind)) /norm(Ca.horGx(Ca.h.posy(ind), Ca.h.posx(ind))));
        end
        if Ca.horGy(Ca.h.posy(ind), Ca.h.posx(ind)) == 0
            Ca.h.posy(ind+1) = Ca.h.posy(ind);
        else
            Ca.h.posy(ind+1) = Ca.h.posy(ind) + round( Ca.horGy(Ca.h.posy(ind), Ca.h.posx(ind)) / norm(Ca.horGy(Ca.h.posy(ind), Ca.h.posx(ind))));
        end
        
        if ind > 1 && (Ca.h.posx(ind)==Ca.h.posx(ind+1) & Ca.h.posy(ind) ==Ca.h.posy(ind+1)) |  (Ca.h.posx(ind-1)==Ca.h.posx(ind+1) & Ca.h.posy(ind) ==Ca.h.posy(ind+1)) | (Ca.h.posx(ind)==Ca.h.posx(ind+1) & Ca.h.posy(ind-1) ==Ca.h.posy(ind+1));
            Ca.h.posx(ind+1) = Ca.h.posx(ind)+randi([-1,1],1);
            Ca.h.posy(ind+1) = Ca.h.posy(ind)+randi([-1,1],1);
        end
        ind = ind +1;
    end
    
    
    Ca.h.closestvesselx = Ca.h.posx(end);
    Ca.h.closestvessely = Ca.h.posy(end);
    
    
    
    %% closest self - distance of both models (vertical, horizontal)    
    for nummodels=1:2
        
        
        Ca.s.posx=[];
        Ca.s.posy=[];
        ind = 1;
        
        if nummodels ==1
            Ca.s.posx(ind)= Ca.v.closestvesselx;
            Ca.s.posy(ind)= Ca.v.closestvessely;
        else
            Ca.s.posx(ind)= Ca.h.closestvesselx;
            Ca.s.posy(ind)= Ca.h.closestvessely;
        end
        
        
        for i=1:Ca.comp.NumObjects
            [I, J] = ind2sub(size(Ca.vert),Ca.comp.PixelIdxList{i});
            if ismembertol(Ca.s.posy,I,5,'DataScale', 1) & ismembertol(Ca.s.posx,J,5,'DataScale',1)
                t.skel = zeros(size(Ca.vert));
                t.ind = 1:size(Ca.comp.PixelIdxList,2);
                t.ind(i)=[];
                for j=1:size(t.ind,2)
                    t.skel(Ca.comp.PixelIdxList{t.ind(j)}) =1;
                end
                break
            end
        end
        
        
        t.D = bwdist(t.skel);
        t.D = imcomplement(t.D);
        [t.Gx ,t.Gy] = imgradientxy(t.D);
        
        e=1;
        while e
            if t.D(Ca.s.posy(ind), Ca.s.posx(ind)) == 1 || t.D(Ca.s.posy(ind), Ca.s.posx(ind)) == 0
                e=0;
                break
            end
            
            if  t.Gx(Ca.s.posy(ind), Ca.s.posx(ind)) == 0
                Ca.s.posx(ind+1) = Ca.s.posx(ind);
            else
                Ca.s.posx(ind+1) = Ca.s.posx(ind) + round( t.Gx(Ca.s.posy(ind), Ca.s.posx(ind)) /norm(t.Gx(Ca.s.posy(ind), Ca.s.posx(ind))));
            end
            if t.Gy(Ca.s.posy(ind), Ca.s.posx(ind)) == 0
                Ca.s.posy(ind+1) = Ca.s.posy(ind);
            else
                Ca.s.posy(ind+1) = Ca.s.posy(ind) + round( t.Gy(Ca.s.posy(ind), Ca.s.posx(ind)) / norm(t.Gy(Ca.s.posy(ind), Ca.s.posx(ind))));
            end
            
            
            if ind > 1 && (Ca.s.posx(ind)==Ca.s.posx(ind+1) & Ca.s.posy(ind) ==Ca.s.posy(ind+1)) |  (Ca.s.posx(ind-1)==Ca.s.posx(ind+1) & Ca.s.posy(ind) ==Ca.s.posy(ind+1)) | (Ca.s.posx(ind)==Ca.s.posx(ind+1) & Ca.s.posy(ind-1) ==Ca.s.posy(ind+1))
                Ca.s.posx(ind+1) = Ca.s.posx(ind)+randi([-1,1],1);
                Ca.s.posy(ind+1) = Ca.s.posy(ind)+randi([-1,1],1);
            end
            ind = ind +1;
            
        end
        selfD(nummodels) = length(Ca.s.posx);
    end
    
    
    
    %% closest vessel in Histo
    % vertical histo model  
    H.v.posx=[];
    H.v.posy=[];
    
    
    ind = 1;
    H.v.posx(ind)= Ca.v.closestvesselx;
    H.v.posy(ind)= Ca.v.closestvessely;
    
    e=1;
    while e
        if H.vertD(H.v.posy(ind), H.v.posx(ind)) == 1 || H.vertD(H.v.posy(ind), H.v.posx(ind)) == 0
            e=0;
            break
        end
        
        if  H.vertGx(H.v.posy(ind), H.v.posx(ind)) == 0
            H.v.posx(ind+1) = H.v.posx(ind);
        else
            H.v.posx(ind+1) = H.v.posx(ind) + round( H.vertGx(H.v.posy(ind), H.v.posx(ind)) /norm(H.vertGx(H.v.posy(ind), H.v.posx(ind))));
        end
        if H.vertGy(H.v.posy(ind), H.v.posx(ind)) == 0
            H.v.posy(ind+1) = H.v.posy(ind);
        else
            H.v.posy(ind+1) = H.v.posy(ind) + round( H.vertGy(H.v.posy(ind), H.v.posx(ind)) / norm(H.vertGy(H.v.posy(ind), H.v.posx(ind))));
        end
        if ind > 1 && (H.v.posx(ind)==H.v.posx(ind+1) & H.v.posy(ind) ==H.v.posy(ind+1)) |  (H.v.posx(ind-1)==H.v.posx(ind+1) & H.v.posy(ind) ==H.v.posy(ind+1)) | (H.v.posx(ind)==H.v.posx(ind+1) & H.v.posy(ind-1) ==H.v.posy(ind+1));
            H.v.posx(ind+1) = H.v.posx(ind)+randi([-1,1],1);
            H.v.posy(ind+1) = H.v.posy(ind)+randi([-1,1],1);
        end
        ind = ind +1;
    end
    
    if length(H.v.posx) >= selfD(1)
        H.v.closestvesselx = nan;
        H.v.closestvessely = nan;
    else
        H.v.closestvesselx = H.v.posx(end);
        H.v.closestvessely = H.v.posy(end);
    end
    
    
    % horizontal histo model
    H.h.posx=[];
    H.h.posy=[];
    
    ind = 1;
    H.h.posx(ind)= Ca.h.closestvesselx;
    H.h.posy(ind)= Ca.h.closestvessely;
    
    ind = 1;
    e=1;
    while e
        if H.horD(H.h.posy(ind), H.h.posx(ind)) == 1 || H.horD(H.h.posy(ind), H.h.posx(ind)) == 0
            e=0;
            break
        end
        
        if  H.horGx(H.h.posy(ind), H.h.posx(ind)) == 0
            H.h.posx(ind+1) = H.h.posx(ind);
        else
            H.h.posx(ind+1) = H.h.posx(ind) + round( H.horGx(H.h.posy(ind), H.h.posx(ind)) /norm(H.horGx(H.h.posy(ind), H.h.posx(ind))));
        end
        if H.horGy(H.h.posy(ind), H.h.posx(ind)) == 0
            H.h.posy(ind+1) = H.h.posy(ind);
        else
            H.h.posy(ind+1) = H.h.posy(ind) + round( H.horGy(H.h.posy(ind), H.h.posx(ind)) / norm(H.horGy(H.h.posy(ind), H.h.posx(ind))));
        end
        
        if ind > 1 && (H.h.posx(ind)==H.h.posx(ind+1) & H.h.posy(ind) ==H.h.posy(ind+1)) |  (H.h.posx(ind-1)==H.h.posx(ind+1) & H.h.posy(ind) ==H.h.posy(ind+1)) | (H.h.posx(ind)==H.h.posx(ind+1) & H.h.posy(ind-1) ==H.h.posy(ind+1))
            H.h.posx(ind+1) = H.h.posx(ind)+randi([-1,1],1);
            H.h.posy(ind+1) = H.h.posy(ind)+randi([-1,1],1);
        end
        ind = ind +1;
        
    end
    
    
    if length(H.h.posx) >= selfD(2)
        H.h.closestvesselx = nan;
        H.h.closestvessely = nan;
    else
        H.h.closestvesselx = H.h.posx(end);
        H.h.closestvessely = H.h.posy(end);
    end
    
    
    %% find segmented elements in images and choose bigger ones        
    for i=1:H.v.comp.NumObjects
        [I, J] = ind2sub(size(Histo.vesselimgbwskel),H.v.comp.PixelIdxList{i});
        if ismembertol(H.v.closestvesselx,J,1,'DataScale', 1) & ismembertol(H.v.closestvessely,I,1,'DataScale', 1)
            H.v.Cind = i;
            break
        end
    end
    
    
    for i=1:Ca.v.comp.NumObjects
        [I, J] = ind2sub(size(Histo.vesselimgbwskel),Ca.v.comp.PixelIdxList{i});
        if ismembertol(Ca.v.closestvesselx,J,1,'DataScale', 1) & ismembertol(Ca.v.closestvessely,I,1,'DataScale', 1)
            Ca.v.Cind = i;
            break
        end
    end
    
    for i=1:H.h.comp.NumObjects
        [I, J] = ind2sub(size(Histo.vesselimgbwskel),H.h.comp.PixelIdxList{i});
        if ismembertol(H.h.closestvesselx,J,1,'DataScale', 1) & ismembertol(H.h.closestvessely,I,1,'DataScale', 1)
            H.h.Cind = i;
            break
        end
    end
    
    
    for i=1:Ca.h.comp.NumObjects
        [I, J] = ind2sub(size(Histo.vesselimgbwskel),Ca.h.comp.PixelIdxList{i});
        if ismembertol(Ca.h.closestvesselx,J,1,'DataScale', 1) & ismembertol(Ca.h.closestvessely,I,1,'DataScale', 1)
            Ca.h.Cind = i;
            break
        end
    end
  
    
    D_v = sqrt((H.v.closestvesselx-Ca.v.closestvesselx)^2 + (H.v.closestvessely-Ca.v.closestvessely)^2);
    D_h = sqrt((H.h.closestvesselx-Ca.h.closestvesselx)^2 + (H.h.closestvessely-Ca.h.closestvessely)^2);
    
    if D_v < D_h
        
        X=imfuse(Ca.vert,H.vert);
        F = figure('Visible',0);
        imshowpair(X, Ca2cells.soma(:,:,cellind),'blend')
        text(120,40,'Multimodal alignment - offset estimation','Color','w','FontWeight','Bold','FontSize',10);

        hold on
        line([H.v.closestvesselx Ca.v.closestvesselx],[H.v.closestvessely Ca.v.closestvessely], 'Color','red','LineWidth',5)
        
        EucD = D_v;
        vec = [H.v.closestvesselx-Ca.v.closestvesselx,H.v.closestvessely-Ca.v.closestvessely];
    else
        X=imfuse(Ca.hor,H.hor);
        F = figure('Visible',0);
        imshowpair(X, Ca2cells.soma(:,:,cellind),'blend')
        text(120,40,'Multimodal alignment - offset estimation','Color','w','FontWeight','Bold','FontSize',10);
        hold on
        line([H.h.closestvesselx Ca.h.closestvesselx],[H.h.closestvessely Ca.h.closestvessely], 'Color','red','LineWidth',5)
        EucD = D_h;
        vec = [H.h.closestvesselx-Ca.h.closestvesselx ,H.h.closestvessely-Ca.h.closestvessely];
    end
    
    Multireg.vesseloffset.EucD(cellind) = EucD;
    Multireg.vesseloffset.vec(cellind,:) = vec;
    t = getframe(gcf);
    close(F);    
    M(:,:,:,cellind)=t.cdata;
    cellind=cellind+1;
    close
end
% Show image with the vessel offset distances
Mmax = max(M,[],4);
figure, imshow(Mmax)
%%
Theta = zeros(1,size(Multireg.vesseloffset.vec,1));
Rho = zeros(1,size(Multireg.vesseloffset.vec,1));
for i=1:size(Multireg.vesseloffset.vec,1)
    [Theta(i), Rho(i)] = cart2pol(Multireg.vesseloffset.vec(i,1),Multireg.vesseloffset.vec(i,2));
end

%% local offsets

%% calculate local offset of cells based on blood vessels

% get min max distance of cells from Ca2+ imaging
minmaxCa2D = max(nanmin(Ca2cells.somaeucD,[],2));

% get Ca2+ soma centroids and try to cluster them
X =[Ca2cells.somacentroid(:,1,:),Ca2cells.somacentroid(:,2,:)];
X=squeeze(X);
X=X';
for clusterrep = 1:size(X,1)
    
    [idx,c] = kmeans(X,clusterrep);
   
    for i=1:max(idx)
        D{i} = euclideanD([c(i,1) c(i,2)],[X(find(idx==i),1) X(find(idx==i),2)]);
    end
    
    Dist=[];
    for i=1:size(D,2)
        Dist =  [Dist; D{i}];
    end
    
    if max(Dist) <= minmaxCa2D
        sprintf(['Number of Clusters: ' num2str(clusterrep) ])
        break
    end
end
Multireg.cluster.assignmentids = idx;
figure, imshow(max(Ca2cells.soma ,[],3))
hold on
text(140,40,'In vivo Ca2+ - Cell cluster detection','Color','w','FontWeight','Bold','FontSize',10);
markers = {'x','.','O','+','*','s','d'};
colors = {'black','red','blue','cyan','green','magenta'};

L = cell(1,max(idx));
for i=1:max(idx)
    
    colind1 = mod(i,6);
    colind2 = mod(i,5);
    if colind1 == 0
        colind1 = 1;
    end
    if colind2 ==0
        colind2 = 2;
    end
    plot(X(idx==i,1),X(idx==i,2),markers{colind1},'color',colors{colind2});
    L{i} = ['Cluster ' num2str(i)];
end
legend(L)

for i=1:max(idx)
    % cartesian coordinate vectors
    currentx =  Multireg.vesseloffset.vec(find(idx==i),1);
    currenty =  Multireg.vesseloffset.vec(find(idx==i),2);
    Multireg.vesseloffset.len(i) = norm(Multireg.vesseloffset.vec(i,:));
    Multireg.cluster.clusteroffsetvec(i,:) = [ nanmean(currentx) nanmean(currenty)];
    Multireg.cluster.clusteroffsetlength(i) = norm([ nanmean(currentx) nanmean(currenty)]);
end

%% find possible matching partners
[Histocell,Calcium_cells, Cellmatch]= multimodal_cell_match_pairs(stack(:,:,1:vesselimgdepth), Histo.segmentedcells(:,:,1:vesselimgdepth), Ca2cells, floor(Ca2cells.Meancellradius));
%% find best matching partner
[Cellmatch] = multimodal_cell_match(Cellmatch, Ca2cells.Meancellradius, Multireg.cluster.clusteroffsetlength, Multireg.cluster.assignmentids);
Calcium_cells.meancellradius = Ca2cells.Meancellradius;

output.cellmatchindex.Ca2 = find(Cellmatch.cellmatchindex);
output.cellmatchindex.histo = Cellmatch.cellmatchindex(find(Cellmatch.cellmatchindex));

output.Calcium_cells = Calcium_cells;
output.Histocell = Histocell;
output.Cellmatch = Cellmatch;
end