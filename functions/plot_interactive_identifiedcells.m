function plot_interactive_identifiedcells(Histology,Cellmatch, Multireg)
% (imagestack,histoperimeter,cellimages, calciumperimeter, matchinindex)
% This function shows automatically matched in vivo recorded cells in post
% hoc histology.
%
% Inputs:

%           Histology : Data structure containing histology imaging data.
%
% Usage :
%           Arrow up : show previous image in histology stack
%      
%           Arrow down: show next image in histology stack
%
%           + : jump + 10 images in histology stack
%
%           - : jump -10 images in histology stack
%
%           Arrow left : show next cell recorded in in vivo Ca2+ imaging
%
%           Arrow right : show previous cell recorded in in vivo Ca2+ imaging
%
%           y : toggle perimeter of in vivo recorded cell on / off
%
%           x : toggle perimeter of post hoc matched cell on / off
%
%           c : show Ca2+ image of in vivo recorded cell
%
%           v : toggle perimeters of all in vivo recorded cells on / off
%
%           b : re-create figure in case user interactio nslows down
%
%           1 - 9 : show histology images of corresponding channel number 
%
%           q : quit
%
%           ESC : quit
%
% Function is written by Philip Anner (2020)

%imagestack = Histologystack;

% Determine if rotation optimization was performed
if Multireg.viewpoint.values.rotx ~= 0 && Multireg.viewpoint.values.roty ~= 0
    if isfield(Histology,'rot')
        histoperimeter = Histology.rot.segmentedcells;
        if isfield(Histology.rot,'Ch')
            imagestacks = Histology.rot.Ch;
        elseif isfield(Histology,'ch')
            imagestacks = Histology.rot.ch;
        else
            error('No histology image stacks found! Please check the data structure of histology imaging data!')
        end
    else
        error('Rotation optimization used but no blood vessels in rotated histology structure found! Run createrotatedstack(Histology, Multireg) first')
    end
else
    histoperimeter = Histology.segmentedcells;
    if isfield(Histology,'Ch')
        imagestacks = Histology.Ch;
    elseif isfield(Histology,'ch')
        imagestacks = Histology.ch;
    else
        error('No histology image stacks found! Please check the data structure of histology imaging data!')
    end
end


histologychannel = 1;
maxhistologychannels = numel(imagestacks);
cellimages = Cellmatch.Calcium_cells.images;
calciumperimeter = Cellmatch.Calcium_cells.soma;
ICind = Cellmatch.cellmatchindex.Ca2;
histocellind = Cellmatch.cellmatchindex.histo;
tmpstack = zeros(size(histoperimeter));
currenthistocell = [];
f1=figure;

% convert to double 
if ~isa(cellimages,'double')
    cellimages = im2double(cellimages);
end

for i = 1:numel(imagestacks)
    if ~isa(imagestacks{i}.stack,'double')
        imagestacks{i}.stack = im2double(imagestacks{i}.stack);
    end
end

e = 1;
togglehistoperim=1;
togglecaperim=1;
toggleimage=0;
toggleallperim=0;
i=1;
caI=1;
hI=1;
change=1;

m = msgbox({'User interaction: ' , 'Arrow up : show previous image in histology stack' ,...
'Arrow down: show next image in histology stack', '+ : jump + 10 images in histology stack',...
'- : jump -10 images in histology stack','Arrow left : show next cell recorded in in vivo Ca2+ imaging',...
'Arrow right : show previous cell recorded in in vivo Ca2+ imaging', 'y : toggle perimeter of in vivo recorded cell on / off',...
'x : toggle perimeter of post hoc matched cell on / off', 'c : show Ca2+ image of in vivo recorded cell',...
'v : toggle perimeters of all in vivo recorded cells on / off', 'b : re-create figure in case user interactio nslows down',...
'1 - 9 : show histology images of corresponding channel number ','q : quit','ESC : quit'});
waitfor(m);

while e
    figure(f1)
    imhandle = imshow(imagestacks{histologychannel}.stack(:,:,i));
    
    [tmpstack, currenthistocell]=plotcurrstackimage(imhandle, imagestacks{histologychannel}.stack,cellimages, histoperimeter,calciumperimeter,i,ICind, caI, histocellind,hI,togglehistoperim, togglecaperim, change,tmpstack,toggleimage, currenthistocell, toggleallperim,histologychannel);
    change=0;
    
    waitforbuttonpress;
    command=gcf;
    val=double(get(command,'CurrentCharacter'));
    try
        switch val
            
            case 49
                % 1 plot histology channel 1
                histologychannel = 1;
            case 50
                % 2 plot histology channel 2
                if 2 <= maxhistologychannels
                    histologychannel = 2;
                end
                
            case 51
                % 3 plot histology channel 3
                if 3 <= maxhistologychannels
                    histologychannel = 3;
                end
            case 52
                % 4 plot histology channel 4
                if 4 <= maxhistologychannels
                    histologychannel = 4;
                end
            case 53
                % 5 plot histology channel 5
                if 5 <= maxhistologychannels
                    histologychannel = 5;
                end
            case 54
                % 6 plot histology channel 6
                if 6 <= maxhistologychannels
                    histologychannel = 6;
                end
                
            case 55
                % 7 plot histology channel 7
                if 7 <= maxhistologychannels
                    histologychannel = 7;
                end
                
            case 56
                % 8 plot histology channel 8
                if 8 <= maxhistologychannels
                    histologychannel = 8;
                end
            case 57
                % 9 plot histology channel 9
                if 9 <= maxhistologychannels
                    histologychannel = 9;
                end
                
           % + key
            case 43
                 i=i+30;
                if i >= size(imagestacks{histologychannel}.stack,3)
                    i=i - 30;
                end
                %- key
            case 45
              i=i-30;
                if i <= 0
                    i=i+30;
                end
                
                
            case 29
                
                caI = caI + 1;
                hI = hI + 1;
                change=1;
                if caI > size(ICind,1)
                    caI = caI-1;
                    hI = hI - 1;
                    change=0;
                end
                
                while hI < size(ICind,1) && histocellind(hI) == 0
                   % caI=caI+1;
                    hI=hI+1;
                end
                  
                % leftArrow
            case 28
                
                caI = caI - 1;
                hI=hI - 1;
                change=1;
                if hI == 0
                    hI = hI+1;
                    caI = caI + 1;
                    change=0;
                end
                
                while hI > 0 && histocellind(hI) == 0
                    %caI=caI-1;
                    hI=hI-1;
                end
                
                %upArrow
            case 30
                i=i-1;
                if i == 0
                    i=i+1;
                end
                
                
                %downArrow
            case 31
                i=i+1;
                if i > size(imagestacks{histologychannel}.stack,3)
                    i=i - 1;
                end
                %case 1
            case 49
                i=1;
                
                
                %case 2
            case 50
                
                i=i-10;
                if i < 1
                    i=i + 10;
                end
                
                %case 3
            case 51
                i=i+10;
                if i > size(imagestacks{histologychannel}.stack,3)
                    i=i - 10;
                end
                
                %case 4
            case 52
                i=size(imagestacks{histologychannel}.stack,3);
                
                % Y - toggle ca perimeter on / off
            case 121
                
                if togglecaperim
                    togglecaperim = 0;
                else
                    togglecaperim = 1;
                end
                
                
                % X - toggle histo perimeter on / off
            case 120
                
                if togglehistoperim
                    togglehistoperim = 0;
                else
                    togglehistoperim = 1;
                end
                
                
                % c - toggle cell image / perimeter
            case 99
                
                if toggleimage
                    toggleimage = 0;
                else
                    toggleimage = 1;
                end
                
                % v - show all ca2+ perim
            case 118
                if toggleallperim
                    toggleallperim = 0;
                else
                    toggleallperim = 1;
                end
                
                % b - new window
            case 98
                close
                f1=figure;
                
                
                % q - Quit
            case 113
                e = 0;
                close
            % ESC - quit    
            case 27    
                e=0;
                close
                
            otherwise
                sprintf('unknown command')
                
        end
    catch
        warning('unknown command');
    end
    
end
end


function [tmpstack, cellrange]=plotcurrstackimage(figurehandle, imagestack,cellimages, histoperimeter,calciumperimeter,i, ICind,caI, histocellind,hI, togglehistoperim, togglecaperim, changeplot,tmpstack,toggleimage, cellrange, toggleallperim,histologychannel)
titlestr=['Imagestack: channel: ' num2str(histologychannel) ', image : ' num2str(i) '/' num2str(size(imagestack,3)) 'Hc: ' num2str(histocellind(hI)) ' Cell: ' num2str(ICind(caI)) ];
set(figurehandle, 'CData', imagestack(:,:,i));

if toggleimage
    set(figurehandle, 'CData', cellimages(:,:,ICind(caI)));
    return
end

if togglehistoperim
    if changeplot
        allregions  = bwconncomp(histoperimeter,26);
        tmpstack = zeros(size(histoperimeter));
        tmpstack(allregions.PixelIdxList{histocellind(hI)})=1;
        pixels=tmpstack;
        
        ind = find(pixels,1,'first');
        [~, ~, cellstart ]=ind2sub(size(pixels),ind);
        
        ind = find(pixels,1,'last');
        [~, ~, cellend]=ind2sub(size(pixels),ind);
        
        
        cellrange = ['. Cell range: ' num2str(cellstart) ' - ' num2str(cellend)];
    else      
        X = gca;
        visboundaries(X,(tmpstack(:,:,i)),'Color','b');
    end
end

if togglecaperim
    X = gca;
    visboundaries(X,(calciumperimeter(:,:,ICind(caI))),'Color','r');
end

if toggleallperim
    s=size(calciumperimeter);
    perimeter=zeros(s(1:2));
    perimeterindex = zeros(size(calciumperimeter));
    for i=1:size(calciumperimeter,3)
        
        if ~isa(calciumperimeter,'logical')
            calciumperimeter(:,:,i) = imbinarize(calciumperimeter(:,:,i));
        end
        calciumperimeter(:,:,i)=bwperim(calciumperimeter(:,:,i));
        perimeterindex(:,:,i) = calciumperimeter(:,:,i).*i*5;
    end
    perimeter=max(perimeterindex,[],3);
    perimeter = label2rgb(uint8(perimeter) ,'hsv', [0 0 0]);
    XX = figurehandle.CData;
    AA = imfuse(XX,perimeter ,'blend');
    figurehandle.CData = AA;
end

titlestr = strcat(titlestr, cellrange);
title(titlestr);

end

