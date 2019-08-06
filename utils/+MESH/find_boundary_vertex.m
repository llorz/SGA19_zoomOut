%2018-09-02
% find the boundary vertex id (and boundary edge list)
function [vid,edge_list] = find_boundary_vertex(mesh)
TRIV = mesh.surface.TRIV;
HE = [[TRIV(:,1),TRIV(:,2)];...
    [TRIV(:,2),TRIV(:,3)];...
    [TRIV(:,3),TRIV(:,1)]]; % half-edges

E = [HE; HE(:,[2,1])];
[val,~,ic] = unique(E,'rows');
a_counts = accumarray(ic,1);
edge_list = val(a_counts == 1,:);
edge_list = edge_list(ismember(edge_list,HE,'rows'),:);
vid = unique(reshape(edge_list,[],1));
end