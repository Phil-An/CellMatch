function [output_txt, G]= GraphCursorCallback(obj,event_obj,val,p,r1,r2,fig)
% This callback function enables interactive modifications to the graph
% structure containing tracked cells across experiments.
%
%
% Function is written by Philip Anner (2020)

% initialize output variable
G = getappdata(fig,'Graph');
% Get selected Node
h = get(event_obj,'Target');
pos = get(event_obj,'Position');
ind = find(h.XData == pos(1) & h.YData == pos(2), 1);
output_txt = {['Node ' num2str(ind)]};

switch val
    
    case 1
        % show image of a cell
        implay(imadjust(G.Nodes.Images{ind}, [0 0.6], [0 1 ]));
        title([G.Nodes.Name{ind}]);
    case 2
        % show cell soma perimeter
        figure, imshow(G.Nodes.Perimeter{ind});
        title([G.Nodes.Name{ind}]);
    case 3
        % show total cell and soma cell perimeter
        figure
        imshowpair(G.Nodes.Perimeter{ind},bwperim(G.Nodes.Area{ind}));
        title([G.Nodes.Name{ind}]);
    case 4
        % Remove edge
        G.Nodes.Name{ind};
        idx=(strfind(G.Edges.EndNodes,G.Nodes.Name{ind}));
        [r, ~ ]= find(not(cellfun('isempty', idx)));
        for i=1:size(r,1)
            str{i}=  strcat(G.Edges.EndNodes{r(i),1},' - ', G.Edges.EndNodes{r(i),2});
        end
        
        figure
        edge_handle = uicontrol('Style', 'popup',...
            'String', str,...
            'Position', [20 340 100 50]);
        uicontrol('Style', 'pushbutton', 'String', 'Delete Edge','Callback', {@clearedge,edge_handle,r,r1,r2,fig});
    case 5
        % Add edge
        str = G.Nodes.Name;
        figure
        edge_handle = uicontrol('Style', 'popup',...
            'String', str,...
            'Position', [20 340 100 50]);
        uicontrol('Style', 'pushbutton', 'String', 'Add Edge','Callback', {@addedgegraph,edge_handle,ind,r1,fig});
    case 6
        % Delete perimeter
        G.Nodes.Perimeter{ind}=ones(size(G.Nodes.Perimeter{1}));
        G.Nodes.Area{ind}=ones(size(G.Nodes.Perimeter{1}));
        G.Nodes.Soma{ind}=ones(size(G.Nodes.Perimeter{1}));
        setappdata(fig,'Graph',G);
        
    case 7
        % Delete Node
        G = rmnode(G, ind);
        setappdata(fig,'Graph',G);
        if r1.Value == 1
            p=plot(G,'ArrowSize',15,'NodeLabel',G.Nodes.Name, 'MarkerSize', 10);
            layout(p,'subspace','UseGravity',true,'Iterations',0);
        else
            p=plot(G,'ArrowSize',15,'NodeLabel',1:size(G.Nodes,1) ,'MarkerSize', 10);
            layout(p,'subspace','UseGravity',true,'Iterations',0);
        end
        wbc = centrality(G,'indegree');
        wbco = centrality(G,'outdegree');
        tmp = wbc./wbco;
        tmp(find(isnan(tmp)))=1;
        p.NodeCData= tmp;
        
    case 8
        % plot trace
        figure
        plot(G.Nodes.Traces{ind});
        title([G.Nodes.Name{ind}]);
end

end

% Removing edge
function fig = clearedge(object_handle,event, edge_handle,str,r1,r2,fig)
G = getappdata(fig,'Graph');
del = str(edge_handle.Value);
G=rmedge(G,del);
setappdata(fig,'Graph',G);
close

if r1.Value == 1
    p=plot(G,'ArrowSize',15,'NodeLabel',G.Nodes.Name, 'MarkerSize', 10);
else
    p=plot(G,'ArrowSize',15,'NodeLabel',1:size(G.Nodes,1) ,'MarkerSize', 10);
end
wbc = centrality(G,'indegree');
wbco = centrality(G,'outdegree');
tmp = wbc./wbco;
tmp(find(isnan(tmp)))=1;
p.NodeCData= tmp;
end

%Adding edge
function  fig = addedgegraph(object_handle,event, edge_handle,r1,r2,fig);
G = getappdata(fig,'Graph');
G = addedge(G,r1,edge_handle.Value,1);
setappdata(fig,'Graph',G);
close

if r2.Value == 1
    p=plot(G,'ArrowSize',15,'NodeLabel',G.Nodes.Name, 'MarkerSize', 10);
else
    p=plot(G,'ArrowSize',15,'NodeLabel',1:size(G.Nodes,1) ,'MarkerSize', 10);
end
wbc = centrality(G,'indegree');
wbco = centrality(G,'outdegree');
tmp = wbc./wbco;
tmp(find(isnan(tmp)))=1;
p.NodeCData= tmp;
end