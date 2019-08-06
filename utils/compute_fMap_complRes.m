% use the complex resolvent mask to compute a functional map
function [C12, Dlb] = compute_fMap_complRes(S1,S2,B1,B2,Ev1,Ev2,fct_src,fct_tar,para,mask_type)
numEigsSrc = size(B1,2); numEigsTar = size(B2,2);
%--------------------------------------------------------------------------
% Descriptors
assert(size(fct_src,2)==size(fct_tar,2));
% Normalization
no = sqrt(diag(fct_src'*S1.A*fct_src))';
fct_src = fct_src ./ repmat(no, [S1.nv,1]);
no = sqrt(diag(fct_tar'*S2.A*fct_tar))';
fct_tar = fct_tar ./ repmat(no, [S2.nv,1]);
%--------------------------------------------------------------------------
% Multiplication Operators
numFct = size(fct_src,2);
OpSrc = cell(numFct,1);
OpTar = cell(numFct,1);
for i = 1:numFct
    OpSrc{i} = B1'*S1.A*(repmat(fct_src(:,i), [1,numEigsSrc]).*B1);
    OpTar{i} = B2'*S2.A*(repmat(fct_tar(:,i), [1,numEigsTar]).*B2);
end
Fct_src = B1'*S1.A*fct_src;
Fct_tar = B2'*S2.A*fct_tar;
% Orientation-preserving Operators
compute_all_OrientationOp = @(S,B,fct) ...
    cellfun(@(f) OrientationOp(S,B,f),mat2cell(fct,size(fct,1),ones(size(fct,2),1)),'un',0);
F11_all = compute_all_OrientationOp(S1,B1,fct_src);
F22_all = compute_all_OrientationOp(S2,B2,fct_tar);
%% all energy terms and the corresponding gradient
max_ev = max(max(Ev1),max(Ev2));
Ev1 = Ev1./max_ev;
Ev2 = Ev2./max_ev;
switch mask_type
    case 'standard' % the standard Laplacian Mask
        Dlb = (repmat(Ev2, [1,numEigsSrc]) - repmat(Ev1', [numEigsTar,1])).^2;
        Dlb = Dlb/norm(Dlb, 'fro');
    case 'complRes' % the complex resolvent Laplacian Mask
        sigma = 0.5;
        Ev1 = Ev1.^sigma;
        Ev2 = Ev2.^sigma;
        
        % complex part
        x1 = Ev2./(1 + Ev2.^2);
        x2 = Ev1'./(1+ Ev1'.^2);
        Dlb1 = (repmat(x1, [1,numEigsSrc]) - repmat(x2, [numEigsTar,1])).^2;
        
        % real part
        x1 = 1./(1 + Ev2.^2);
        x2 = 1./(1+ Ev1'.^2);
        Dlb2 = (repmat(x1, [1,numEigsSrc]) - repmat(x2, [numEigsTar,1])).^2;
        % equal weight
        Dlb = Dlb1 + Dlb2;
        Dlb = Dlb/norm(Dlb,'fro');
    case 'slant'   % the slanted mask presented in "Partial functional map"
        sigma = 0.03;
        est_rank = sum(Ev1-max(Ev2)<0);
        Dlb = zeros(numEigsTar,numEigsSrc);
        for i=1:numEigsTar
            for j=1:numEigsSrc
                slope = est_rank/numEigsTar;
                direction = [1 slope];
                direction = direction./norm(direction);
                Dlb(i,j) = exp(-sigma*sqrt(i.^2 + j.^2))*norm(cross([direction 0], [i,j, 0]-[1 1 0]));
            end
        end
        Dlb = Dlb/norm(Dlb, 'fro');
end

% orientation term
funOrient_direct = @(C) sum(cellfun(@(F11,F22) 0.5*norm(C*F11 - F22*C, 'fro')^2,F11_all, F22_all));
gradOrient_direct = @(C) sum(cell2mat(cellfun(@(F11,F22) reshape(C*(F11*F11') - F22'*C*F11 - F22*C*F11' + F22'*F22*C,[],1),...
    F11_all, F22_all,'un',0)),2);
% descriptors term
funDesp = @(C) 0.5*norm(C*Fct_src - Fct_tar,'fro')^2;
gradDesp = @(C) reshape((C*Fct_src - Fct_tar)*Fct_src',[],1);
% commutativity with descriptors
funCommDesp = @(C) sum(cell2mat(cellfun(@(X,Y) 0.5*norm(X*C - C*Y,'fro')^2, OpTar', OpSrc', 'UniformOutput', false)), 2);
gradCommDesp = @(C) sum(cell2mat(cellfun(@(X,Y) reshape(X'*(X*C - C*Y) - (X*C - C*Y)*Y',[],1), OpTar', OpSrc', 'UniformOutput', false)), 2);
% commutativity with LB
funCommLB = @(C) sum(sum((C.^2 .* Dlb)/2));
gradCommLB =@(C) reshape(C.*Dlb,[],1);

%% Fmap Computation
%--------------------------------------------------------------------------
constFct = sign(B1(1,1)*B2(1,1))*[sqrt(sum(S2.area)/sum(S1.area)); zeros(numEigsTar-1,1)];
F_lb = zeros(numEigsTar*numEigsSrc, 1); F_lb(1) = constFct(1);

C = mat_projection(eye(numEigsTar,numEigsSrc)); % identity initialization

Dlb1 = (repmat(Ev2, [1,numEigsSrc]) - repmat(Ev1', [numEigsTar,1])).^2;
Dlb1 = Dlb1/norm(Dlb1, 'fro'); % make the weight consistent across different masks
% make the relative weight consistent w.r.t different terms
a = para.a/funDesp(C); b = para.b/funCommDesp(C); alpha = para.alpha/funOrient_direct(C);
c = para.c/(sum(sum((C.^2 .* Dlb1)/2)));

myObj = @(C) a*funDesp(C) + b*funCommDesp(C) + c*funCommLB(C) + alpha*funOrient_direct(C);
myGrad = @(C) a*gradDesp(C) + b*gradCommDesp(C) + c*gradCommLB(C) + alpha*gradOrient_direct(C);
funObj = @(F) deal(myObj(reshape(F,numEigsTar, numEigsSrc)),...
    myGrad(reshape(F,numEigsTar, numEigsSrc)));


funProj = @(F) [constFct; F(numEigsTar+1:end)];
options.maxIter = 500;
options.verbose = 1;
tic;
C12 = reshape(minConf_PQN(funObj, F_lb, funProj, options), [numEigsTar,numEigsSrc]);
t = toc; fprintf('done %.4fs.\n', t);
end


function [C,v,eval] = mat_projection(W)
n1 = size(W,1);
n2 = size(W,2);
[s,v,d] = svd(full(W));
C = s*eye(n1,n2)*d';
v = diag(v);
eval = [range(v),mean(v),var(v)];
end