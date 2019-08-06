function desc_projected = descriptors_projection(desc, B, A)
% project the given descriptors into the given basis
% desc: a set of descriptors defined on S
% B: the basis of S where the fMap lives in
% A: the area matrix of the mesh S
if nargin > 2
    B_inv = B'*A;
else
    B_inv = pinv(B);
end
desc_projected = B_inv * desc;
end