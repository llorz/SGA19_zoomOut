% by Adrien Poulenard :D
function [ SymOp, HOp ] = OrientationOp( surface, Lb , H)
%F Summary of this function goes here
%   Detailed explanation goes here

% compute heat kernel signature
t = 100; % numTimes
if nargin < 3 % no input descriptors
    H = hks( surface, t );
end
% compute its gradient
G = face_grads(surface, H);

% normalize it so that it has unit norm
vn  = sqrt(sum(G'.^2))';
vn = vn + eps;
vn = repmat(vn,1,3);
G = G./vn;

% rotate it by pi/2 using cross product with the normal
JGsrc = -J(surface,G);

% create 1st order differential operators associated with the vector fields
SymOp = Lb'*surface.A*vf2op(surface, JGsrc)*Lb;
HOp = Lb'*surface.A*vf2op(surface, G)*Lb;
end

%%
function [ h ] = hks( shape, t )
%HKS Summary of this function goes here
%   Detailed explanation goes here

LB = shape.evecs(:, 1:200 );
eigs = shape.evals(1:200 );
squaredLb = LB.^2;
T = exp(-abs(eigs*t));

h = squaredLb*T;
end
%%
function G = face_grads(mesh, f)

X = mesh.surface.VERT;
T = mesh.surface.TRIV;
Nf = mesh.normals_face;

v1 = X(T(:,1),:);
v2 = X(T(:,2),:);
v3 = X(T(:,3),:); 

Ar = 0.5*sum((cross(v3-v2,v1-v2)).^2,2); % area per face

Ar = repmat(Ar, 1, 3);
G = repmat(f(T(:,1)),1,3).*cross(Nf, v3 - v2)./(2*Ar) + ...
    repmat(f(T(:,2)),1,3).*cross(Nf, v1 - v3)./(2*Ar) + ...
    repmat(f(T(:,3)),1,3).*cross(Nf, v2 - v1)./(2*Ar);
end
%%
function [ JVf ] = J( shape, Vf )
%J Summary of this function goes here
%   Detailed explanation goes here
JVf = cross(Vf, shape.normals_face, 2);
end
%%
function op = vf2op(mesh, Vf)

X = mesh.surface.VERT;
T = mesh.surface.TRIV;
nf = mesh.nf;
nv = mesh.nv;

Nf = mesh.normals_face;

v1 = X(T(:,1),:);
v2 = X(T(:,2),:);
v3 = X(T(:,3),:); 
C1 = v3 - v2;
C2 = v1 - v3;
C3 = v2 - v1;
Jc1 = cross(Nf, C1);
Jc2 = cross(Nf, C2);
Jc3 = cross(Nf, C3);

I = [T(:,1);T(:,2);T(:,3)];
J = [T(:,2);T(:,3);T(:,1)];
Sij = 1/6*[dot(Jc2,Vf,2); dot(Jc3,Vf,2); dot(Jc1,Vf,2)];
Sji = 1/6*[dot(Jc1,Vf,2); dot(Jc2,Vf,2); dot(Jc3,Vf,2)];
In = [I;J;I;J];
Jn = [J;I;I;J];
Sn = [Sij;Sji;-Sij;-Sji];
W = sparse(In,Jn,Sn,nv,nv);

M = mass_matrix(mesh);
%op = spdiags(1./mesh.origAreaWeights,0,nv,nv)*W;
op = spdiags(1./sum(M,2),0,nv,nv)*W;
end
%%
function M = mass_matrix(mesh)

T = mesh.surface.TRIV; 
X = mesh.surface.VERT;
v1 = X(T(:,1),:);
v2 = X(T(:,2),:);
v3 = X(T(:,3),:); 

Ar = 0.5*sum((cross(v3-v2,v1-v2)).^2,2); % area per face
nv = mesh.nv;

I = [T(:,1);T(:,2);T(:,3)];
J = [T(:,2);T(:,3);T(:,1)];
Mij = 1/12*[Ar; Ar; Ar];
Mji = Mij;
Mii = 1/6*[Ar; Ar; Ar];
In = [I;J;I];
Jn = [J;I;I];
Mn = [Mij;Mji;Mii];
M = sparse(In,Jn,Mn,nv,nv);
end
