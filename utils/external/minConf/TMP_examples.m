clear all
close all

%% Non-negative Least Squares
% Solve min_w (1/2)norm(X*w-y)^2, s.t. w >= 0.

nInstances = 1000;
nVars = 100;
X = randn(nInstances,nVars);
w = randn(nVars,1);
y = X*w + randn(nInstances,1);

funObj = @(w)SquaredError(w,X,y);
LB = zeros(nVars,1);
UB = inf(nVars,1);

fprintf('Solving non-negative least-squares problem...\n');
w = minConf_TMP(funObj,zeros(nVars,1),LB,UB);
stem(w);title('Non-negative least-squares coefficients');
pause

%% L1-Regularization Least Squares
% Solve min_w (1/2)norm(X*w-y)^2 + lambda*sum(abs(w))
%  by formulating as
%       min_w (1/2)norm([X -X]*w - y)^2 + lambda*sum(w), s.t. w >= 0

lambda = 1000*ones(nVars,1);

regObj = @(w)nonNegGrad(w,lambda,funObj);
LB = zeros(2*nVars,1);
UB = inf(2*nVars,1);

fprintf('Solving L1-regularized least-squares problem...\n');
w = minConf_TMP(regObj,zeros(2*nVars,1),LB,UB);
w = w(1:nVars)-w(nVars+1:end);
stem(w);title('L1-regularized least-squares coefficients');
pause

%% Logistic Regression with Bounded Coefficients
% Solve min_w sum_i log(1 + exp(-y(i)*w*X(i,:)), s.t. -1 <= w <= 1

y = sign(y);

funObj = @(w)LogisticLoss(w,X,y);
LB = -ones(nVars,1);
UB = ones(nVars,1);

fprintf('Solving bounded logistic regression problem...\n');
w = minConf_TMP(funObj,zeros(nVars,1),LB,UB);
stem(w);title('Bounded logistic regression coefficients');
pause


%% Dual Support Vector Machines (no bias or regularized bias)
% Solve min_alpha (1/2)alpha'*A*alpha - sum(alpha), s.t. 0 <= alpha <= C,
%   where A_ij = y_iy_jK(x_i,x_j), and we use the RBF kernel with sigma = 1.

% Generate data
nInstances = 250;
[X,y] = makeData('classificationNonlinear',nInstances,2,2);

sigma = 1;
K = kernelRBF(X,X,sigma);
A = diag(y)*K*diag(y);

C = 1;
funObj = @(alpha)dualSVMLoss_noBias(alpha,A,y);
LB = zeros(nInstances,1);
UB = C*ones(nInstances,1);

fprintf('Solving dual SVM problem (no bias)...\n');
alpha = minConf_TMP(funObj,zeros(nInstances,1),LB,UB);
stem(alpha);title('Dual SVM coefficients');

yhat = sign(sum((diag(alpha.*y)*K)))';
trainErr = sum(yhat~=y)/numel(y)
pause

%% Ordinal Logistic Regression
% Solve min_{w,gamma} -sum(log(F(gamma(y+1) - X*w) - F(gamma(y) - X*w))),
%  where F(x) = 1/(1+exp(-x), -inf < 0 < gamma(1) < gamma(2) < ... < inf

% Generate Data
nInstances = 1000;
nVars = 10;
nClasses = 5;
X = randn(nInstances,nVars);
w = randn(nVars,1);
gamma = sort(randn(nClasses-1,1));
z = X*w;
y = zeros(nInstances,1);
y(z < gamma(1)) = 1;
for class = 2:nClasses-1
    y(z >= gamma(class-1) & z < gamma(class)) = class;
end
y(z >= gamma(nClasses-1)) = nClasses;

% Standardize columns and add bias
X = standardizeCols(X);
X = [ones(nInstances,1) X];
nVars = nVars+1;

% First try Multinomial Logistic
fprintf('Training multinomial logistic classifier\n');
model = classificationSoftmax(X,y,struct('nClasses',nClasses));
yhat = model.predictFunc(model,X);
trainErr_MLR = sum(abs(yhat-y));

% Ordinal Logistic
funObj = @(w)OrdinalLogisticLoss2(w,X,y,nClasses);
LB = [-inf(nVars,1);zeros(nClasses-2,1)];
UB = inf(nVars+nClasses-2,1);
w_init = zeros(nVars,1);
gamma_init = sort(rand(nClasses-2,1));
fprintf('Training ordinal logistic classifier\n');
wGamma = minConf_TMP(funObj,[w_init;gamma_init],LB,UB);
w = wGamma(1:nVars);
gamma = [-inf;0;cumsum(wGamma(nVars+1:end));inf];

% Predict labels on training data
z = X*w;
yhat = zeros(nInstances,1);
for c = 1:nClasses
   yhat(z > gamma(c)) = c;
end
trainErr_OLR = sum(abs(yhat-y));

fprintf('Training error of multinomial logistic regression: %f\n',trainErr_MLR/nInstances);
fprintf('Training error of ordinal logistic regression: %f\n',trainErr_OLR/nInstances);
pause

%% Kernel Ordinal Logistic Regression

% Generate data
nInstances = 500;
nVars = 2;
nClasses = 5;
X = randn(nInstances,nVars);
nExamplePoints = 3;
examplePoints = randn(nExamplePoints,nVars);
thresholds = [0;cumsum(2*rand(nClasses-1,1))];
y = zeros(nInstances,1);
for i = 1:nInstances
    dists = sum((repmat(X(i,:),nExamplePoints,1) - examplePoints).^2,2);
    y(i,1) = max(find(min(dists) > thresholds));
end
X = [ones(nInstances,1) standardizeCols(X)];

nTrain = nInstances/2;
Xtrain = X(1:nTrain,:);
ytrain = y(1:nTrain);
Xtest = X(nTrain+1:end,:);
ytest = y(nTrain+1:end);

% First try kernel multinomial logistic
sigma = 1;
lambda = 1e-5;
model = classificationKernelSoftmax(Xtrain,ytrain,struct('nClasses',nClasses,'kernelFunc',@kernelRBF,'kernelArgs',sigma,'lambda',lambda));
yhat = model.predictFunc(model,Xtest);
testErr_KMLR = sum(yhat~=ytest)/length(ytest);
testDist_KMLR = sum(abs(yhat-ytest))/length(ytest);

% Now try kernel ordinal logistic
Ktrain = kernelRBF(Xtrain,Xtrain,sigma);

% Set up problem
w = zeros(nTrain,1);
gamma = ones(nClasses-2,1);
LB = [-inf(nTrain,1);zeros(nClasses-2,1)];
UB = inf(nTrain+nClasses-2,1);
funObj_sub = @(w)OrdinalLogisticLoss2(w,Ktrain,ytrain,nClasses);
funObj = @(w)penalizedKernelL2_subset(w,Ktrain,1:nTrain,funObj_sub,lambda);

% Solve optiamization
wGamma = minConf_TMP(funObj,[w;gamma],LB,UB);
w = wGamma(1:nTrain);
gamma = [-inf;0;cumsum(wGamma(nTrain+1:end));inf];

% Predict on test data
Ktest = kernelRBF(Xtest,Xtrain,sigma);
z = Ktest*w;
yhat = zeros(size(ytest));
for c = 1:nClasses
   yhat(z > gamma(c)) = c; 
end
testErr_KOLR = sum(yhat~=ytest)/length(ytest);
testDist_KOLR = sum(abs(yhat-ytest))/length(ytest);

fprintf('Test error of kernel multinomial logistic regression: %f (distance = %f)\n',testErr_KMLR,testDist_KMLR);
fprintf('Test error of kernel ordinal logistic regression: %f (distance = %f)\n',testErr_KOLR,testDist_KOLR);
pause

%% Graphical LASSO
% Solve min_{W positive-definite} logdet(sigma+W), s.t. |W_ij| <= lambda

% Load data
load 20news_w100.mat
docs = full(documents)';
[nSamples,nVars] = size(docs);
mu = mean(docs);
centered = docs-repmat(mu,16242,1);
sigma = (1/nSamples)*(centered'*centered);

% Set up and solve problem
lambda = .01;
funObj = @(W)logdetFunction(W,sigma);
LB = -lambda*ones(nVars);
UB = lambda*ones(nVars);
W = eye(nVars);
W(:) = minConf_TMP(funObj,W(:),LB(:),UB(:));
K = inv(sigma+W);
K(abs(W) < lambda) = 0;
clf;
drawGraph(K~=0,'labels',wordlist);
pause

%% Associative Conditional Random Fields (trained with pseudo-likelihood)
% Optimize pseudo-likelihood in CRF with Ising potentials, subject to
% constraint that edges are sub-modular (and hence the optimal MAP can be
% found using graph cuts)

% Load Data
load X.mat
y = int32(1+X);
X = X + randn(size(X))/2;
[nRows,nCols] = size(X);
nNodes = nRows*nCols;
nStates = 2;
y = reshape(y,[1 1 nNodes]);
X = reshape(X,1,1,nNodes);

% Set up problem in UGM
adj = latticeAdjMatrix(nRows,nCols);
edgeStruct = UGM_makeEdgeStruct(adj,nStates);
tied = 1;
Xnode = [ones(1,1,nNodes) UGM_standardizeCols(X,tied)];
sharedFeatures = [1 0];
Xedge = UGM_makeEdgeFeaturesInvAbsDif(Xnode,edgeStruct.edgeEnds,sharedFeatures);
ising = 1;
[nodeMap,edgeMap,w] = UGM_makeCRFmaps(Xnode,Xedge,edgeStruct,ising,tied);
nParams = length(w);
funObj = @(w)UGM_CRF_PseudoNLL(w,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct);
UB = inf(nParams,1);
LB = [-inf;-inf;0;0];
w = minConf_TMP(funObj,w,LB,UB)

[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
MAP = UGM_Decode_GraphCut(nodePot,edgePot,edgeStruct);
imagesc(reshape(X,nRows,nCols));colormap gray
figure
imagesc(reshape(MAP,nRows,nCols));colormap gray
