% original: by Anne Gehre
% Implemenation of "Interactive Curve Constrained Functional Maps", by
% Anne Gehre, Michael Bronstein, Leif Kobbelt, and Justin Solomon (SGP
% 2018)

function [vertexids,dist] = find_shortest_path(S,v1,v2,type)

if nargin < 4, type = 'Euclidean'; end

mesh = S.connectivity;
mesh.X = S.surface.VERT;
switch type
    case 'Euclidean'
        vertexids = shortestPathEuklidean(mesh, v1, v2);
    otherwise
        error(['Unsupported type: ',type]);
end

vertexids = reshape(vertexids,[],1);
X = S.surface.VERT;
edge = [vertexids(1:end-1), vertexids(2:end)];
dist = sum(sqrt(sum((X(edge(:,1),:) - X(edge(:,2),:)).^2,2)));
end


function vertexids = shortestPathEuklidean(mesh, id0, id1)
X1 = mesh.X(mesh.E2V(:,1),:);
X2 = mesh.X(mesh.E2V(:,2),:);
dist = X1 - X2;
norms = zeros(size(dist, 1), 1);
for i=1:size(dist,1)
    norms(i) = norm(dist(i,:));
end

G = graph(mesh.E2V(:,1)',mesh.E2V(:,2)', norms');

[vertexids, ~] = shortestpath(G,id0, id1);
end

%% compute_curvature.m is missing: check with Anne
% for Cmax..
function vertexids = shortestPathCurvature(mesh, id0, id1)
X1 = abs(mesh.Cmax(mesh.E2V(:,1),:));
X2 = abs(mesh.Cmax(mesh.E2V(:,2),:));
curvatureM = (X1 + X2);
maxC = max(curvatureM);
curvature = repmat(maxC, size(curvatureM));

X1 = mesh.X(mesh.E2V(:,1),:);
X2 = mesh.X(mesh.E2V(:,2),:);
dist = X1 - X2;
norm = repmat((sum(dist.^2, 2)).^0.5, 1, 3);
ndist = dist./norm;
Cdir = mesh.Umin(mesh.E2V(:,1),:);
z = sum(ndist.*(Cdir),2);
curvatureDev = acosd(z);
indices = find(curvatureDev > 90);
curvatureDev(indices) = 180-curvatureDev(indices);
G = graph(mesh.E2V(:,1)',mesh.E2V(:,2)', (curvatureDev').*(curvature'));

[vertexids, ~] = shortestpath(G,id0, id1);
end

%% uncomment to use anisotropic distances (function shortestCurvatureTensor.m)
% [mesh.Umin,mesh.Umax,mesh.Cmin,mesh.Cmax,mesh.Cmean,mesh.Cgauss,mesh.Normal] = compute_curvature(mesh.X,mesh.T);
% mesh.Umin = mesh.Umin';
% mesh.Umax = mesh.Umax';
% mesh.Normal = mesh.Normal';
%mesh.ShapeOperator = shapeOperator(mesh);%
%mesh.anisotropicDistance = anisotropicDistance(mesh);

%
% function ad = anisotropicDistance(mesh)
%     X1 = mesh.X(mesh.E2V(:,1),:);
%     X2 = mesh.X(mesh.E2V(:,2),:);
%     dist = X1 - X2;
%     for i=1:length(mesh.ShapeOperator)
%         ad(i) = (abs(dist(i,:)*(mesh.ShapeOperator{i}*dist(i,:)'))).^0.5;
%     end
% end

%%
% function vertexids = shortestPathGeodesic(mesh, id0, id1)
% options.method='discrete';
% [D,S,Q] = perform_fast_marching_mesh(mesh.X, mesh.T, id0, options);
% [path, vertexids, plist] = compute_geodesic_mesh(D,mesh.X, mesh.T, id1, options);
% end