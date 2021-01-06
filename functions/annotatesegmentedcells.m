function segmentedcells = annotatesegmentedcells(imagestack,Histostack_EQ, varargin)
% This function let's a user manually annotate cells in histology volumes.
% Manual user annotation is computer assisted. Clicking in a region of a
% cell is usually sufficient to automatically detect the shape of the
% neuron.
%
% Input : 
%           imagestack : a 3D matrix of a laser confocal microscopy 
%                        histology scan of brain tissue. 
%
%           Histostack_EQ : a histogram-equalized and pre-processed version
%                           of the histology laser confocal microscopy scans
%
%           VARARGIN(1) : if you have already started annotating cells, you
%                         can define previously segmented cells as an input
%                         variable and modify them:
%                         segmentedcells = annotatesegmentedcells(imagestack,Histostack_EQ, segmentedcells)
%
% Output:
%           segmentedcells : 3D matrix of segmented neurons in the histology.
%
% Functionality :
%
%   [Arrow RIGHT] : show next image in histology volume
%
%   [Arrow LEFT] : show previous image in histology volume
%   
%   [Arrow UP] : jump -10 images in histology volume
%
%   [Arrow DOWN] : jump +10 images in histology volume
%
%   [Space] + mouse left click on a cell: detect and add the 3D structure
%                of a cell
%
%   [Space] + mouse right button click on a cell: delete a 2D shape of a
%                segmented cell (on a single image of the histology volume)
%
%   [Space] + mouse middle button click on a cell: delete a 3D shape of a
%                segmented cell (on multiple images of the histology volume)
%
%   [y] : toggle segmented cell perimeters on / off
%
%   [x] : rebuild figure in case the user interaction slows down
%
%   [c] +  mouse left click on a cell: detect and add a 2D structure
%                of a cell
%
%   [c] + mouse right button click on a cell: detect and add a small 2D
%                structure of a cell (more restrictive than [c] +  mouse left click)
%
%   [c] + mouse middle button click on a cell: detect and add a very small 2D
%                structure of a cell (more restrictive than [c] +  mouse right click)
%
%   [v] : freehand drawing tool for annotating cell
%
%   [b] : undo last action
%
%   [a] : show segmented cell perimeters from the previous image of the
%         histology stack
%
%   [w] : show segmented cell perimeters from the following image of the
%         histology stack
%
%   [a] : separate touching annotated cells by draginw a line
%
%
% Function is written by Philip Anner (2020)

 m= msgbox({'Functionality : ', ...
'[Arrow RIGHT] : show next image in histology volume', ...
'[Arrow LEFT] : show previous image in histology volume', ...
'[Arrow UP] : jump -10 images in histology volume', ...
'[Arrow DOWN] : jump +10 images in histology volume', ...
'[Space] + mouse left click on a cell: detect and add the 3D structure of a cell', ...
'[Space] + mouse right button click on a cell: delete a 2D shape of a segmented cell (on a single image of the histology volume)', ...
'[Space] + mouse middle button click on a cell: delete a 3D shape of a segmented cell (on multiple images of the histology volume)', ...
'[y] : toggle segmented cell perimeters on / off', ...
'[x] : rebuild figure in case the user interaction slows down', ...
'[c] +  mouse left click on a cell: detect and add a 2D structure of a cell', ...
'[c] + mouse right button click on a cell: detect and add a small 2D structure of a cell (more restrictive than [c] +  mouse left click)', ...
'[c] + mouse middle button click on a cell: detect and add a very small 2D structure of a cell (more restrictive than [c] +  mouse right click)', ...
'[v] : freehand drawing tool for annotating cell', ...
'[b] : undo last action', ...
'[a] : show segmented cell perimeters from the previous image of the histology stack', ...
'[w] : show segmented cell perimeters from the following image of the histology stack', ...
'[a] : separate touching annotated cells by draginw a line'});

waitfor(m);

if numel(varargin) == 1
    segmentedcells = varargin{1};
elseif numel(varargin) > 2
   error('Wrong number of input arguments! Please call annotatesegmentedcells(imagestack,Histostack_EQ) or annotatesegmentedcells(imagestack,Histostack_EQ, segmentedcells)') 
else
    segmentedcells = imbinarize(zeros(size(imagestack)));
end

if ~isa(segmentedcells,'logical')
    segmentedcells = imbinarize(segmentedcells);
end

tic;
numberdeleted=0;
numberadded=0;
previousperim=0;
nextperim=0;
cellind = 1;
Histostack_EQ=im2uint8(mat2gray(Histostack_EQ));
initialregions= regionprops(segmentedcells,'PixelList');
initialcells = size(initialregions,1);

f1=figure;
i=1;
e=1;
showperim = 1;


while e
    
    figure(f1)
    imhandle = imshow(imagestack(:,:,i));
    plotcurrstackimage(imagestack, segmentedcells,i,showperim,previousperim,cellind, nextperim);
    xlabel('SPACE -left click: add 3D -right click: delete 2D -middle button: delete 3D, Y: perim on/off, X: new figure, B: undo, A: previous perim, C -> left: ad 2D cell -right: add small 2D -middle: add very small 2D')
    
    waitforbuttonpress;
    command=gcf;
    val=double(get(command,'CurrentCharacter'));
    
    %viscircles to prevent figure resize
    viscircles(([0 0 ]), 1,'EdgeColor','b');
    try
        switch val
            
            case 28
                %right Arrow
                i=i-1;
                if i == 0
                    i=i+1;
                end
                
            case 29
                %left Arrow
                i=i+1;
                if i > size(Histostack_EQ,3)
                    i=i - 1;
                end
                
            case 30
                %up Arrow
                i=i-10;
                if i <= 0
                    i=i+10;
                end

            case 31
                %down Arrow
                i=i+10;
                if i > size(Histostack_EQ,3)
                    i=i - 10;
                end
                
            case 32
                % space - select 1 element
                tI = segmentedcells(:,:,i);
                tI = (tI==cellind);
                regions = regionprops(tI,'PixelList');
                [x,y,button] = ginput(1);

                if button ==1
                    % left button add cell
                    [~, J]=regionGrowing(Histostack_EQ,[round(y) round(x) i],((max(max(Histostack_EQ(:,:,i))))- min(min(Histostack_EQ(:,:,i))))*0.1, 10);%, true,'tfSimpl');
                    previoussegmentedcells=segmentedcells;
                    J= J.*cellind;
                    segmentedcells=max(segmentedcells,J);
                    numberadded=numberadded+1;
                end
                
                if button == 2
                    %middle mouse button - delete 3D object
                    tV = (segmentedcells==cellind);
                    allregions = regionprops(tV,'PixelList');
                    segmentedcellindex=1;
                    foundcell=0;
                    for j=1:size(allregions,1)
                        if ismember([round(x) round(y) i],allregions(j).PixelList,'rows' )
                            foundcell=1;
                            break
                        end
                        segmentedcellindex=segmentedcellindex+1;
                    end
                    
                    if foundcell == 1
                        previoussegmentedcells = tV;
                        imgind = allregions(segmentedcellindex).PixelList;
                        
                        for j=1:size(imgind,1)
                            tV(imgind(j,2), imgind(j,1),imgind(j,3))=0;
                        end
                        segmentedcells(segmentedcells==cellind)=tV;
                        numberdeleted=numberdeleted+1;
                    end
                end
                
                if button == 3
                    %right button delete 2D cell
                    segmentedcellindex=1;
                    foundcell=0;
                    for j=1:size(regions,1)
                        if ismember([round(x) round(y)],regions(j).PixelList,'rows' )
                            foundcell=1;
                            break
                        end
                        segmentedcellindex=segmentedcellindex+1;
                    end
                    
                    if foundcell == 1
                        previoussegmentedcells=segmentedcells;
                        imgind = regions(segmentedcellindex).PixelList;
                        tI(imgind(:,2), imgind(:,1))=0;
                        tI = segmentedcells(:,:,i).*tI;
                        
                        tI2 = (segmentedcells(:,:,i) > 0 & segmentedcells(:,:,i)~= cellind);
                        tI2 = segmentedcells(:,:,i).*tI2;
                        segmentedcells(:,:,i) = max(tI2, tI);
                        numberdeleted=numberdeleted+1;
                    else
                        sprintf('No cell found!')
                    end
                end

            case 121
                % Y - toggle perimeter on / off
                if showperim
                    showperim = 0;
                else
                    showperim = 1;
                end
                
            case 120
                % X - new figure
                close
                f1 = figure;
                
            case 98
                % b - for undo
                segmentedcells=previoussegmentedcells;
                
            case 97
                % a - for previous perimeter
                if previousperim
                    previousperim = 0;
                else
                    previousperim = 1;
                end
                % s for applying previous perimeters
            case 115
                if nextperim
                    nextperim = 0;
                else
                    nextperim = 1;
                end

            case 119
                % w - for next image perimeter overlap
                previoussegmentedcells=segmentedcells;
                if i >1
                    segmentedcells(:,:,i) = max(segmentedcells(:,:,i),segmentedcells(:,:,i-1));
                end
                
                
            case 101
                % s for applying next image perimeters
                previoussegmentedcells=segmentedcells;
                if i < size(Histostack_EQ,3)
                    segmentedcells(:,:,i) = max(segmentedcells(:,:,i),segmentedcells(:,:,i+1));
                end
                
            case 99
                % c - region growing on 2D image
                sprintf('2D image - region growing')
                
                tI = segmentedcells(:,:,i);
                tI = (tI==cellind);
                regions = regionprops(tI,'PixelList');
                
                [x,y,button] = ginput(1);
                % left button add cell
                if button ==1
                    [~, J]=regionGrowing(Histostack_EQ(:,:,i),[round(y) round(x)],((max(max(Histostack_EQ(:,:,i))))- min(min(Histostack_EQ(:,:,i))))*0.1, 10);%, true,'tfSimpl');
                    J=J.*cellind;
                    % store previous matrix for undo
                    previoussegmentedcells=segmentedcells;
                    % calculate new 2d image
                    tI=max(tI,J);
                    % update segmentatio n stack
                    segmentedcells(:,:,i)=max(segmentedcells(:,:,i),tI);
                    %counter
                    numberadded=numberadded+1;
                end
                
                % right button - add small cell area
                if button == 3
                    [~, J]=regionGrowing(Histostack_EQ(:,:,i),[round(y) round(x)],((max(max(Histostack_EQ(:,:,i))))- min(min(Histostack_EQ(:,:,i))))*0.05, 10);%, true,'tfSimpl');
                    J=J.*cellind;
                    % store previous matrix for undo
                    previoussegmentedcells=segmentedcells;
                    % calculate new 2d image
                    tI=max(tI,J);
                    % update segmentation stack
                    segmentedcells(:,:,i)=max(segmentedcells(:,:,i),tI);
                    %counter
                    numberadded=numberadded+1;  
                end
                
                % middle button - add very small cell area
                if button == 2
                    [~, J]=regionGrowing(Histostack_EQ(:,:,i),[round(y) round(x)],((max(max(Histostack_EQ(:,:,i))))- min(min(Histostack_EQ(:,:,i))))*0.01, 10);%, true,'tfSimpl');
                    J=J.*cellind;
                    % store previous matrix for undo
                    previoussegmentedcells=segmentedcells;
                    % calculate new 2d image
                    tI=max(tI,J);
                    % update segmentation stack
                    segmentedcells(:,:,i)=max(segmentedcells(:,:,i),tI);
                    %counter
                    numberadded=numberadded+1;
                end
                
                
            case 118
                % v -- draw free hand roi
                hold off
                imshow(imagestack(:,:,i))
                
                tI = segmentedcells(:,:,i);
                tI = (tI==cellind);
                % store previous matrix for undo
                previoussegmentedcells=segmentedcells;
                
                % draw new cell
                h = imfreehand
                J = h.createMask.*cellind;
                tI=max(tI,J);
                
                % update segmentation stack
                segmentedcells(:,:,i)=max(segmentedcells(:,:,i),tI);
                %counter
                numberadded=numberadded+1;
                hold on

            case 100
                % s - separate touching cells with line tool
                hold off
                imshow(imagestack(:,:,i))
                
                tI = segmentedcells(:,:,i);
                tI = (tI==cellind);
                tI_allother = (segmentedcells(:,:,i) ~= cellind & segmentedcells(:,:,i) > 0);
                tI_allother = tI_allother .* 9;
                
                % store previous matrix for undo
                previoussegmentedcells=segmentedcells;
                
                % draw line and mask image
                h = imline;
                J = h.createMask;
                J=double(J);
                J = imgaussfilt(double(J),1);
                J = imbinarize(J);
                tI =tI.*~J;
                tI = max(tI, tI_allother);

                % update segmentation stack
                segmentedcells(:,:,i)=tI;
                hold on
 
            case 113
                % q - Quit
                e = 0;
                close
                printsummary(initialcells, numberadded,numberdeleted,segmentedcells);
                
                
            case 27
                % ESC - Quit
                e = 0;
                close
                 printsummary(initialcells, numberadded,numberdeleted,segmentedcells);
                
            otherwise
                sprintf('unknown command')
        end
        
    catch
        warning('unknown command');
    end
    
end
if ~islogical(segmentedcells)
    segmentedcells = imbinarize(segmentedcells);
end
end


function plotcurrstackimage( Histostack_EQ, segmentedcells,i, toggleperim,previousperim,cellind, nextperim)
title(['Image ' int2str(i) ' of ' num2str(size(Histostack_EQ,3)) ' CI: ' num2str(cellind)])
if toggleperim
    hold on
    perimeter = label2rgb(bwperim(segmentedcells(:,:,i)) .*segmentedcells(:,:,i),'hsv', [0 0 0]);
    p=rgb2gray(perimeter);
    p=(p>1);
    h=imshow(perimeter);
    set(h, 'AlphaData', p)
end

if (previousperim && i>1)
    hold on
    visboundaries((segmentedcells(:,:,i-1)),'Color','r');
end

if (nextperim && i<size(Histostack_EQ,3))
    hold on
    visboundaries((segmentedcells(:,:,i+1)),'Color','g');
end

end

function printsummary(initialcells, numberadded,numberdeleted,segmentedcells)
timespent=toc;

disp('====SUMMARY====');
disp(['Initial cells: ' num2str(initialcells)]);
disp(['Cells added: ' num2str(numberadded)]);
disp(['Cells deleted: ' num2str(numberdeleted)]);

if ~islogical(segmentedcells)
    segmentedcells = imbinarize(segmentedcells);
end
finalregions= regionprops(segmentedcells,'PixelList');
finalcells = size(finalregions,1);
disp(['Cells finally: ' num2str(finalcells)]);
disp(['Time spent: ' num2str(timespent/60) ' minutes']);
end