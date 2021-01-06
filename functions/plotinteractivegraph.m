function Graph = plotinteractivegraph(G)
% This function allows a user to interactively assess and modify cells
% tracked across experiments. Segmented cell shapes can be removed for
% individual cells.
%
% Input:
%           G : Graph structure, derived from construct_graph() function.
%
% Function is written by Philip Anner (2020)

% initialize output variable
Graph = [];
Graph = G;
% plot Graph structure
f = figure('CloseRequestFcn',@closecallback);
p= plot(Graph,'ArrowSize',15,'NodeLabel',Graph.Nodes.Name, 'MarkerSize', 10);
setappdata(f,'Graph',Graph);

% colorize relation of in-/ out edges
wbc = centrality(Graph,'indegree');
wbco = centrality(Graph,'outdegree');
tmp = wbc./wbco;
tmp(find(isnan(tmp)))=1;
p.NodeCData= tmp;

% initiate interactive plotting
bg = uibuttongroup('Parent', f,'Visible','off',...
    'Position',[0 0 .2 1],...
    'SelectionChangedFcn',@namecallback);
r1 = uicontrol('Parent',bg,'Style',...
    'radiobutton',...
    'String','Names',...
    'Position',[10 350 100 30],...
    'HandleVisibility','off');

r2 = uicontrol('Parent',bg,'Style','radiobutton',...
    'String','IDs', 'Position',[10 250 100 30], 'HandleVisibility','off');
bg.Visible = 'on';
popup = uicontrol('Style', 'popup',...
    'String', {'Image','Perimeter - soma','Perimeter - all cell',...
    'Delete Edges','Add Edge','Delete Perimeter','Delete Node','Show Trace'},...
    'Position', [100 600 100 50],'Callback', @graphcallback);

hdt = datacursormode;
hdt.UpdateFcn = @(obj,event_obj) GraphCursorCallback(obj,event_obj,popup.Value,p, r1,r2,f);

    % Callback to switch between cell names and IDs 
    function namecallback(object_handle, event)
        if strcmp(event.NewValue.String, 'Names')
            p=plot(Graph,'ArrowSize',15,'NodeLabel',Graph.Nodes.Name, 'MarkerSize', 10);
        else
            p=plot(Graph,'ArrowSize',15,'NodeLabel',1:size(Graph.Nodes,1) ,'MarkerSize', 10);
        end
        
        wbc = centrality(Graph,'indegree');
        wbco = centrality(Graph,'outdegree');
        tmp = wbc./wbco;
        tmp(find(isnan(tmp)))=1;
        p.NodeCData= tmp;
    end

    % Callback for closing and exporting modified graph structure
    function closecallback(src,callbackdata)
        selection = questdlg('Done?',...
            'Closing cell tracking validation',...
            'Yes','No','Yes');
        switch selection
            case 'Yes'
                Graph = getappdata(callbackdata.Source,'Graph');
                delete(src)
            case 'No'
                return
        end
    end

waitfor(f)
end