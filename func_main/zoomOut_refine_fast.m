function T12_refined = zoomOut_refine_fast(S1, S2, T12, para, seed)

k = para.k_init;

if nargin > 4, rng(seed); end

X1_all = S1.evecs(euclidean_fps(S1.surface, para.num_samples), :);
X2_all = S2.evecs(euclidean_fps(S2.surface, para.num_samples), :);

B1 = S1.evecs(:,1:k);
B2 = S2.evecs(:,1:k);
C21 = B1\B2(T12, :);
map12 = knnsearch(X2_all(:, 1:k)*C21', X1_all(:, 1:k));

[~, C21] = zoomOut_refine(X1_all, X2_all, map12, para);

T12_refined = knnsearch(S2.evecs(:,1:size(C21,2))*C21', S1.evecs(:,1:size(C21,1)));
end


% Euclidean farthest point sampling.
function idx = euclidean_fps(surface,k,seed)

C = [surface.X surface.Y surface.Z];
nv = size(C,1);

if(nargin<3)
    idx = randi(nv,1);
else
    idx = seed;
end

dists = bsxfun(@minus,C,C(idx,:));
dists = sum(dists.^2,2);

for i = 1:k
    maxi = find(dists == max(dists));
    maxi = maxi(1);
    idx = [idx; maxi];
    newdists = bsxfun(@minus,C,C(maxi,:));
    newdists = sum(newdists.^2,2);
    dists = min(dists,newdists);
end

idx = idx(2:end);
end