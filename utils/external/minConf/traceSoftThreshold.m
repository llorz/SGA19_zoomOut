function [w] = traceSoftThreshold(w,alpha,lambda,n,p)
W = reshape(w,n,p);
[U,S,V] = svd(W);
W = U*max(0,S-lambda*alpha)*V';
w = W(:);