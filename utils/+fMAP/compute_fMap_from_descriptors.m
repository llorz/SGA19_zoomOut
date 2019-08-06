function C12 = compute_fMap_from_descriptors(S1, S2, fct1, fct2, para)
import fMAP.*
if nargin < 5
    para.type_orient = 'direct';
    para.a = 2e-1; para.b = 1e-2;
    para.c = 8e-4; para.d = 1e-1;
    para.fMap_size = [50, 50];
end

k1 = para.fMap_size(1);
k2 = para.fMap_size(2);

B1 = S1.evecs(:, 1:k1); Ev1 = S1.evals(1:k1);
B2 = S2.evecs(:, 1:k2); Ev2 = S2.evals(1:k2);


% normalize and project the corresponding descriptors
fct1 = descriptors_normalization(fct1, S1.A);
fct2 = descriptors_normalization(fct2, S2.A);
Fct1 = descriptors_projection(fct1, B1, S1.A);
Fct2 = descriptors_projection(fct2, B2, S2.A);

% compute the descriptor operators
compute_all_DescOp = @(B, A, fct)...
    cellfun(@(f) B'*A*(repmat(f, [1, size(B,2)]).*B),...
    mat2cell(fct,size(fct,1),ones(size(fct,2),1)),'un',0);

DescOp1 = compute_all_DescOp(B1, S1.A, fct1);
DescOp2 = compute_all_DescOp(B2, S2.A, fct2);

% compute the orientation operators
compute_all_OrientOp = @(S,B,fct) ...
    cellfun(@(f) OrientationOp(S,B,f),mat2cell(fct,size(fct,1),ones(size(fct,2),1)),'un',0);
OrientOp1 = compute_all_OrientOp(S1, B1, fct1);
OrientOp2 = compute_all_OrientOp(S2, B2, fct2);


% --------------------------------------------------------------------------
%                    Collect all the regularizers
% --------------------------------------------------------------------------
[fval_desc, grad_desc] = ...
    regularizer_descriptor_preservation(Fct1, Fct2);

switch para.type_orient
    case 'direct'
        [fval_comm_orientOp, grad_comm_orientOp] = ...
            regularizer_operator_commutativity(OrientOp1, OrientOp2);
    case 'symmetric'
        [fval_comm_orientOp, grad_comm_orientOp] = ...
            regularizer_operator_commutativity(OrientOp1, OrientOp2, true);
end

[fval_comm_descOp, grad_comm_descOp] = ...
    regularizer_operator_commutativity(DescOp1, DescOp2);

[fval_comm_LB, grad_comm_LB] = ...
    regularizer_laplacian_commutativity(Ev1, Ev2);

% --------------------------------------------------------------------------
%                       Fmap Computation
% --------------------------------------------------------------------------
% the first col of the fMap is fixed by the first eigenvectors
constFct = sign(B1(1,1)*B2(1,1))*[sqrt(sum(S2.area)/sum(S1.area)); zeros(k2-1,1)];
F_lb = zeros(k1*k2, 1); F_lb(1) = constFct(1);

C = mat_projection(eye(k2,k1)); % identity initialization
% fix the relative weight of each term w.r.t. the given weights
a = para.a/fval_desc(C); b = para.b/fval_comm_descOp(C);
c = para.c/(fval_comm_LB(C)+1e-4); d = para.d/fval_comm_orientOp(C);

% complete energy and gradient
myObj = @(C12) a*fval_desc(C12) + b*fval_comm_descOp(C12) + c*fval_comm_LB(C12) + d*fval_comm_orientOp(C12);
myGrad = @(C12) a*grad_desc(C12) + b*grad_comm_descOp(C12) + c*grad_comm_LB(C12) + d*grad_comm_orientOp(C12);
funObj = @(C12_vec) deal(myObj(reshape(C12_vec,k2,k1)),...
    myGrad(reshape(C12_vec, k2, k1)));


funProj = @(F) [constFct; F(k2+1:end)];
options.verbose = 1;
fprintf('Optimizing the fMap with %s operator...',para.type_orient);tic;
C12 = reshape(minConf_PQN(funObj, F_lb, funProj, options), [k2, k1]);
t = toc; fprintf('done %.4fs.\n', t);
end