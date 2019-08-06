clear all
close all

%% L1-Regularization
% We first solve min_W ||X*W - Y||_F^2 + lambda*sum_{ij}|W_{ij}|.

n = 500;
p = 100;
k = 50;
X = randn(n,p);
W = diag(rand(p,1) > .9)*randn(p,k) + randn(p,1)*randn(1,k);
Y = X*W + randn(n,k);

% Smooth Part of Objective Function
funObj1 = @(W)SimultaneousSquaredError(W,X,Y);

% Non-Smooth Part
lambda = 500;
funObj2 = @(W)lambda*sum(abs(W));

% Prox Operator
funProx = @(W,alpha)sign(W).*max(abs(W)-lambda*alpha,0);

% Optimize with QNST
fprintf('Optimizing with L1-regularization\n');
W(:) = minConf_QNST(funObj1,funObj2,W(:),funProx);

imagesc(W);
sparsityW = nnz(W)
rowSparsityW = nnz(sum(abs(W),2))
rankW = rank(W)
pause

%% Group L_{1,2}-Regularization
% We now replace the L1-norm regularizer with the group-L1 regularizer,
%   lambda*sum_g ||W_g||_2, where in this case we will define the groups
%   'g' as the rows of the matrix (so each group penalizes an original
%   input feature across the tasks)

% Non-Smooth part
groups = repmat([1:p]',1,k);
nGroups = p;
lambda = 5000*ones(nGroups,1);
funObj2 = @(w)groupL12regularizer(w,lambda,groups);

% Prox operator
funProx = @(w,alpha)groupSoftThreshold(w,alpha,lambda,groups);

% Optimize with QNST
fprintf('Optimizing with Group L1-regularization\n');
W(:) = minConf_QNST(funObj1,funObj2,W(:),funProx);

imagesc(W);
sparsityW = nnz(W)
rowSparsityW = nnz(sum(abs(W),2))
rankW = rank(W)
pause

%% Group L_{1,inf}-Regularization
% We now replacethe L2-norm in the group L1-regularizer by the
% infinity-norm

% Non-Smooth Part
lambda = 20000*ones(nGroups,1);
funObj2 = @(w)groupL1infregularizer(w,lambda,groups);

% Prox operator
funProx = @(w,alpha)groupInfSoftThreshold(w,alpha,lambda,groups);

% Optimize with QNST
fprintf('Optimizing with L1-regularization (inf-norm of groups)\n');
W(:) = minConf_QNST(funObj1,funObj2,W(:),funProx);

imagesc(W);
sparsityW = nnz(W)
rowSparsityW = nnz(sum(abs(W),2))
rankW = rank(W)
pause

%% Combined L1- and Group L1-Regularization
% We now consider applying both the group L1-regularizer to select rows,
% and the regular L1-regularizer to select within the rows

% Non-Smooth Part
lambda1 = 500;
lambda2 = 5000*ones(nGroups,1);
funObj2 = @(w)lambda1*sum(abs(w)) + groupL12regularizer(w,lambda2,groups);

% Prox operator
funProx1 = @(w,alpha,lambda1)sign(w).*max(abs(w)-lambda1*alpha,0);
funProx = @(w,alpha)groupSoftThreshold(funProx1(w,alpha,lambda1),alpha,lambda2,groups);

% Optimize with QNST
fprintf('Optimizing with combined L1- and group L1-regularization\n');
W(:) = minConf_QNST(funObj1,funObj2,W(:),funProx);

imagesc(W);
sparsityW = nnz(W)
rowSparsityW = nnz(sum(abs(W),2))
rankW = rank(W)
pause

%% Nuclear norm-regularization

% Non-Smooth Part
lambda = 1000;
toMat = @(w)reshape(w,p,k);
funObj2 = @(w)lambda*sum(svd(toMat(w)));

% Prox Operator
funProx = @(w,alpha)traceSoftThreshold(w,alpha,lambda,p,k);

% Optimize with QNST
fprintf('Optimizing with Nuclear norm-regularization\n');
W(:) = minConf_QNST(funObj1,funObj2,W(:),funProx);

imagesc(W);
sparsityW = nnz(W)
rowSparsityW = nnz(sum(abs(W),2))
rankW = rank(W)
pause

%% Sparse plus low-rank

% Smooth Part
funObj1 = @(ww)SimultaneousSquaredError(ww,[X X],Y);

% Non-Smooth Part
lambda1 = 500;
lambdaT = 1000;
funObj2 = @(ww)lambda1*sum(abs(ww(1:end/2))) + lambdaT*sum(svd(toMat(ww(end/2+1:end))));

% Prox Operator
fprintf('Optimizing with an L1-regularized matrix plus a nuclear norm-regularized matrix\n');
funProx = @(ww,alpha)[funProx1(ww(1:end/2),alpha,lambda1);traceSoftThreshold(ww(end/2+1:end),alpha,lambdaT,p,k)];
    
% Optimize with QNST
WW = [W W];
WW(:) = minConf_QNST(funObj1,funObj2,WW(:),funProx);
W1 = WW(:,1:k);
W2 = WW(:,k+1:end);

imagesc([W1 W2]);
sparsityWW = [nnz(W1) nnz(W2)]
rowSparsityWW = [nnz(sum(abs(W1),2)) nnz(sum(abs(W2),2))]
rankWW = [rank(W1) rank(W2)]
pause
