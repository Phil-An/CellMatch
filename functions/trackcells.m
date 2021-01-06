function [combinedcells]= trackcells(thresh, varargin)
% This function automatically tracks the same neurons across multiple in
% vivo ca2+ imaging recordings.
%
% Inputs:
%           thresh : Is a threshold for the minimum joint probability density
%                    function of two cells for accepting them as being
%                    tracked across experiments. This parameter depends on
%                    the accurracy of the registration results and the
%                    observed somatic cell size. This parameter must be
%                    estimated empirically.
%           VARARGIN : individual in vivo Ca2+ imaging data structures REC.
%           For aligning multiple recordings, all REC structures must be
%           defined as input parameters.
%
% Outputs:
%           allocation
%           combinedcells
%           probabilities
%           example:
%           [allocation combinedcells probabilities]= trackcells(0.6, REC1, REC2, REC3)
%
%
%
% Function is written by Philip Anner (2020)

sessions=varargin;
%% calculate mean soma radius
nodelcells = 0;
ind=1;

% determine total number of cells recorded
totcells = cellfun(@(x) size(x.IC,3),sessions,'UniformOutput',false);
totcells = cell2mat(totcells);
totcells = sum(totcells);

% preallocate memory for variables
cellsize = zeros(1,totcells);



for i=1:length(sessions)
    for j=1:size(sessions{i}.soma,3)
        
        % check if counter exceeds session size due to deleted cells
        if j + nodelcells > size(sessions{i}.soma,3)
            break
        end
        
        if ~islogical(sessions{i}.soma(:,:,j))
            t= regionprops(imbinarize(sessions{i}.soma(:,:,j)) , 'EquivDiameter');
        else
            t= regionprops(sessions{i}.soma(:,:,j) , 'EquivDiameter');
        end
        
        if size(t,1) == 1
            cellsize(ind) = t.EquivDiameter;
        elseif size(t,1) == 0
            sprintf([ 'No cell found - delete cell: ' num2str(j) ' from session: ' num2str(i)])
            sessions{i}.IC(:,:,j) = [];
            sessions{i}.perim(:,:,j) = [];
            sessions{i}.areas(:,:,j) = [];
            sessions{i}.soma(:,:,j) = [];
            sessions{i}.centroids(:,:,j)= [];
            sessions{i}.names(j) = [];
            sessions{i}.traces.trace(j) = [];
            sessions{i}.traces(j).events = [];
            combinedcells.deletednodes.session(ind) = i;
            combinedcells.deletednodes.node(ind) = j;
            j=j+1;
            ind = ind-1;
            nodelcells = nodelcells +1;
        end
        ind = ind+1;
    end
end
gaussfiltersize = floor(median(cellsize/4));
clear cellsize

%get max length
m=0;
sessionindex=[];
IC = cell(1,length(sessions));
img = cell(1,length(sessions));
orig_length = cell(1,length(sessions));
for i=1:length(sessions)
    IC{i}=squeeze(sessions{i}.centroids);
    IC{i}=round(IC{i})';
    img{i}=sessions{i}.IC;
    orig_length{i}=length(IC{i});
    %get max IC length
    m = max(length(IC{i}),m);
    sessionindex = [ sessionindex; ones(orig_length{i},1).*i];
end



for i=1:(length(sessions))
    names{i}=varargin{i}.names';
end

refsize=size(varargin{1}.IC(:,:,1));
t=zeros((refsize));
t(round(refsize(1)/2), round(refsize(2)/2))=1;
t=imgaussfilt(t,gaussfiltersize); % 10
t=mat2gray(t);
maxPDF=trapz(trapz(t.*t));


for i=1:(length(sessions))
    % bring ICs to same length by padding with 0
    if length(IC{i}) < m
        d = m - length(IC{i});
        IC{i} = [IC{i}; zeros(d,2)];
        img{i} = cat(3,img{i}, zeros([refsize d]));
        pad_diff{i}=d;
    else
        pad_diff{i}=0;
    end
    % create PDF images for m IC sets
    for j=1:m
        t=zeros(refsize);
        
        % if centroid is defined, calculate PDF
        if or(IC{i}(j,1)~=0,IC{i}(j,2)~=0)
            
            
            % check if centroid is out of boundary
            if IC{i}(j,2) > (refsize(1))
                IC{i}(j,2) = (refsize(1));
                
            end
            
            if IC{i}(j,1) > (refsize(2))
                IC{i}(j,1) = (refsize(2));
            end
            
            
            t(IC{i}(j,2), IC{i}(j,1))=1;
            t=imgaussfilt(t,gaussfiltersize);
            PD{i}(:,:,j)=t;
            PD{i}(:,:,j) = mat2gray(PD{i}(:,:,j));
            
        else
            PD{i}(:,:,j)= t;
        end
    end
end

ind=1;
ii=1;
icname = cell(totcells,1);
cellnr = cell(totcells,1);
for i=1:length(names)
    for j=1:size(names{i},2)
        icname{ii,1} =strcat('Session ', num2str(ind), '-',  names{i}(j));
        cellnr{ii,1} = j;
        ii=ii+1;  
    end
    ind=ind+1;
end

res = cell(length(sessions),m);
rescorr = cell(length(sessions),m);
% i..number of IC recordings - reference
for i=1:length(sessions)
    % j .. number of ICs - reference
    for j=1:m
        % k... num of recordings - moving index
        for k=1:length(sessions)
            % k... num of ICs - moving index
            for l=1:m
                % if centroid is defined and ref - moving index is inequal
                if (IC{i}(j,1) ~=0 & IC{i}(j,2) ~=0)  & k ~= i
                    if IC{k}(l,:) ~=[0 0]
                        t =((PD{i}(:,:,j).*PD{k}(:,:,l)));
                        
                        res{i,j}(k,l)=trapz(trapz(t))/maxPDF;
                        
                        t=corr2(img{i}(:,:,j),img{k}(:,:,l));
                        
                        if ~isnan(t)
                            rescorr{i,j}(k,l) =t;
                        else
                            rescorr{i,j}(k,l) = 0;
                        end
                        
                    else
                        res{i,j}(k,l)=0;
                        rescorr{i,j}(k,l) = 0;
                    end
                else
                    res{i,j}(k,l)=0;
                    rescorr{i,j}(k,l) = 0;
                end
            end
        end
    end
end

orig_length_cumulative = cumsum(cell2mat(orig_length));
score = cell(length(sessions),m);
index = 1;
for  i = 1:length(sessions)
    for j=1:m
        [val, ind] =max(res{i,j},[],2);
        [valcorr, indcorr] =max(rescorr{i,j},[],2);
        
        if ~(ind==indcorr)
            dbstop
            sprintf('INCONSISTENT MATCH')
        end
        
        ind(find(val==0))=0;
        score{i,j}.val=valcorr.*(valcorr > thresh);
        score{i,j}.ind=ind.*(val > thresh);
        if size(score{i,j}.ind,1) > 1
            ind = 0;
            for x=1:size(score{i,j}.ind,1)
                if score{i,j}.ind(x)>0
                    if ind > 0
                        indid(x)= score{i,j}.ind(x)+ orig_length_cumulative(ind);
                    else
                        indid(x)= score{i,j}.ind(x);
                    end
                else
                    indid(x)=0;
                end
                
                ind = ind+1;
            end
            score{i,j}.indid=indid;
        end
        index = index+1;
    end
end

index=1;
allocation = cell(totcells, 1);
probabilities = cell(totcells, 1);
allocationid = cell(totcells, 1);
for  i = 1:length(sessions)
    for j=1:m-pad_diff{i}
        allocation{index,1}= score{i,j}.ind;
        probabilities{index,1}= score{i,j}.val;
        allocationid{index,1}= score{i,j}.indid;
        index=index+1;
    end
end

edgestable=[];
weightstable=[];
index=1;
for i=1:length(sessions)
    ref=(index:orig_length_cumulative(i))';

    % get edges
    t=cell2mat(allocationid(index:orig_length_cumulative(i)))';
    t(i,:)=[];
    
    % weights for edges
    tp=cell2mat(probabilities(index:orig_length_cumulative(i)))';
    tp=reshape(tp, [length(sessions),orig_length{i}]);
    tp(i,:)=[];
    
    % IC names
    for j=1:size(t,1)
        inddismiss = find(~t(j,:));
        edges = [ref t(j,:)'];
        edges(inddismiss,:)=[];
        edgestable = [edgestable; edges];
        
        % ADD WEIGHTS HERE
        w = tp(j,:)';
        w(inddismiss,:)=[];
        weightstable = [weightstable; w];
    end
    if i < length(sessions)
        index = orig_length_cumulative(i)+1;
    end
end

% Evaluate if cells were deleted
if ~exist('combinedcells','var')
    combinedcells.deletednodes.session = [];
    combinedcells.deletednodes.node = [];
end

combinedcells.sessionindex = sessionindex;
combinedcells.edges = edgestable;
combinedcells.w= weightstable;
combinedcells.ids= allocationid;
combinedcells.names = icname;
combinedcells.cellnr = cellnr;
end