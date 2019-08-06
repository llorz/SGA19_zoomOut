function desc_normalized = descriptors_normalization(desc, A)
% normalize the descriptor w.r.t. the mesh area
% desc: a set of descriptors defined on mesh S
% A: the area matrix of the mesh S
no = sqrt(diag(desc'*A*desc))';
nv = size(desc,1);
desc_normalized = desc ./ repmat(no, [nv,1]);
end