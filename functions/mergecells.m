function ifenditiedcells = mergecells(H)
% This function computes final representations of cells tracked across
% experiments.
%
% Input:
%           H : Graph structure, derived from construct_graph() function.
%           It is higily recommended to prior manually evaluate and correct
%           matched cells using:
%                               G = plotinteractivegraph(G);
%           followed by calculating the transitive closure by
%                               H = transclosure(G);
%           and removing remove self loops
%                               H = rmedge(H, 1:numnodes(H), 1:numnodes(H));
%           before calling
%                               merge_perim(H)
% Output:
%           ifenditiedcells : structure including paramters of cells
%           tracked across recordings
%
%           ifenditiedcells.allperim : perimeters of total cell areas
%                                      tracked across recordings
%           ifenditiedcells.allarea :  total cell area of cells tracked
%                                      across recordings
%           ifenditiedcells.soma :     cell soma area of cells tracked
%                                      across recordings
%           ifenditiedcells.image :    mean images of cells tracked across
%                                      recordings
%           ifenditiedcells.nmerged :  number indicating the instances of a
%                                      cell matched across recordings
%
%           identifiedcells.recordings.mergedperimeter : movie of total
%           cell perimeters identified in each recording used for matching
%
%           identifiedcells.recordings.mergedsomaperimeter : movie of cell
%           soma perimeters identified in each recording used for matching
%
%
% Function is written by Philip Anner (2020)

% Find all cells and instnaces of a cell across recordings
comp = conncomp(H);
[~,  ia, ~] = unique(comp);

% preallocate memory for variables
allarea = zeros([size(H.Nodes.Area{1}) length(ia)]);
allsoma = zeros([size(H.Nodes.Area{1}) length(ia)]);
allperim = zeros([size(H.Nodes.Area{1}) length(ia)]);
allsomaperim = zeros([size(H.Nodes.Area{1}) length(ia)]);
allmeanimage  = zeros([size(H.Nodes.Area{1}) length(ia)]);
perimlabel = cell(1,length(ia));
nmerged = zeros(1,length(ia));
perimeterindex = zeros([size(H.Nodes.Area{1}) length(ia)]);
somaperimeterindex = zeros([size(H.Nodes.Area{1}) length(ia)]);

index=1;
for i=1:length(ia)
    ind = find(comp==comp(ia(i)));
    
    if size(ind,2) == 1
        % if a cell is identified in 1 recording
        allarea(:,:,index)= logical(H.Nodes.Area{ia(i)});
        allsoma(:,:,index)= logical(H.Nodes.Soma{ia(i)});
        allperim(:,:,index)= logical(bwperim(H.Nodes.Area{ia(i)}));
        allsomaperim(:,:,index)= logical(bwperim(H.Nodes.Soma{ia(i)}));
        allmeanimage(:,:,index)=H.Nodes.Images{ia(i)};
        perimlabel{i}=index;
        nmerged(index)=1;
        index = index+1;
        
    elseif size(ind,2) > 1
        % if a cell is identified in > 1 recordings
        tempimg = zeros(size(H.Nodes.Area{1}));
        for j=1:size(ind,2)
            if j==1
                tempimg = H.Nodes.Area{ind(j)};
                tempimgsoma = H.Nodes.Soma{ind(j)};
                tempimgorig = H.Nodes.Images{ind(j)};
            else
                tempimg = tempimg + H.Nodes.Area{ind(j)} ;
                tempimgsoma = tempimgsoma + H.Nodes.Soma{ind(j)};
                tempimgorig = mean(cat(3,tempimgorig,H.Nodes.Images{ind(j)}),3);
            end
        end
        tempimgsoma= imbinarize(tempimgsoma);
        tempimg= imbinarize(tempimg);
        if ~sum(sum(tempimg))
            sprintf(['Empty IC found! Node: ' num2str(ind(j))])
        end
        allsoma(:,:,index)= tempimgsoma;
        allsomaperim(:,:,index)= bwperim(tempimgsoma);
        allarea(:,:,index)= tempimg;
        allperim(:,:,index)= bwperim(tempimg);
        allmeanimage(:,:,index)= tempimgorig;
        
        perimlabel{i}=ind;
        
        nmerged(index)=j;
        index = index+1;
    end
end

for i=1:size(allperim,3)
    perimeterindex(:,:,i) = allperim(:,:,i).*i;
    somaperimeterindex(:,:,i) = allsomaperim(:,:,i).*i;
end

perimeter=max(perimeterindex,[],3);
mergedperimimage = label2rgb(perimeter ,'hsv', [0 0 0]);

somaperimeter =max(somaperimeterindex,[],3);
mergedsomaperimimage = label2rgb(somaperimeter ,'hsv', [0 0 0]);

ifenditiedcells.recordings.mergedperimeter = mergedperimimage;
ifenditiedcells.recordings.mergedsomaperimeter = mergedsomaperimimage;
ifenditiedcells.allperim=allperim;
ifenditiedcells.allarea=allarea;
ifenditiedcells.soma=allsoma;
ifenditiedcells.nmerged = nmerged;
ifenditiedcells.image = mat2gray(allmeanimage);
end