function G = construct_graph(graph_st, varargin)
% This function creates a structure G to represent cells tracked across
% recordings as a graph.
% Example:
%   G = construct_graph(combinedcells, REC1, REC2, REC3);
%
% Function is written by Philip Anner (2020)

sessiondel = graph_st.deletednodes.session(find(graph_st.deletednodes.session));
nodedel = graph_st.deletednodes.node(find(graph_st.deletednodes.node));

for i=1:size(sessiondel,2)
    SessionIndex= find(graph_st.sessionindex==sessiondel(i));
    deletenode = SessionIndex(nodedel(i));
    graph_st.names(deletenode) = [];
    graph_st.cellnr(deletenode) = [];
end

EdgeTable = table(graph_st.edges, graph_st.w, 'VariableNames',{'EndNodes' 'Weight'});
nnn=cellfun(@(x) char(x), graph_st.names,'UniformOutput', false);

nnn2=cellfun(@(x) double(x), graph_st.cellnr,'UniformOutput', false);
NodeTable = table(nnn, nnn2,graph_st.sessionindex,'VariableNames',{'Name', 'ID', 'Sessionindex'});
clear nnn nnn2
G = digraph(EdgeTable,NodeTable);

IC=cellfun(@(varargin) (squeeze(num2cell(im2double(mat2gray(varargin{:}.IC)),[1 2]))), varargin,'UniformOutput',false );
for i=1:size(sessiondel,2)
    IC{sessiondel(i)}(nodedel(i))=[];
end
G.Nodes.Images = vertcat(IC{:});

Perimeter=cellfun(@(varargin) (squeeze(num2cell(im2double(varargin{:}.perim),[1 2]))), varargin,'UniformOutput',false );
for i=1:size(sessiondel,2)
    Perimeter{sessiondel(i)}(nodedel(i))=[];
end
G.Nodes.Perimeter = vertcat(Perimeter{:});

Area=cellfun(@(varargin) (squeeze(num2cell(im2double(varargin{:}.areas),[1 2]))), varargin,'UniformOutput',false );
for i=1:size(sessiondel,2)
    Area{sessiondel(i)}(nodedel(i))=[];
end
G.Nodes.Area = vertcat(Area{:});

Soma=cellfun(@(varargin) (squeeze(num2cell(im2double(varargin{:}.soma),[1 2]))), varargin,'UniformOutput',false );
for i=1:size(sessiondel,2)
    Soma{sessiondel(i)}(nodedel(i))=[];
end
G.Nodes.Soma = vertcat(Soma{:});

if isfield(varargin{1},'traces')
    Traces = cellfun(@(varargin) (squeeze(num2cell((varargin{:}.traces.trace),[1 2]))), varargin,'UniformOutput',false );
    for i=1:size(sessiondel,2)
        Traces{sessiondel(i)}{1,1}(nodedel(i))=[];
    end
    t=[];
    for i=1:size(Traces,2)
        t = [t;Traces{i}{:}'];
    end
    G.Nodes.Traces = t;
end

G.Edges.EdgeTable=[graph_st.edges graph_st.w];
end