function [func_fval, func_grad] = regularizer_operator_commutativity(Op1, Op2, IfReversing)
% the term to compute the operators via a functional map C12
% Op1 = {Y_1, ..., Y_k} operators defined on S1
% Op2 = {X_1, ..., X_k} corresponding operators defined on S2
% e.g., the operators can be the descriptor commutativity or orientation
% operator
% The regularizer looks like:
%           sum_{i = 1, ..., k} || X_i C12 - C12 Y_i ||_F^2
% when IfReversing is set to true, the regularizer looks like
%           sum_{i = 1, ..., k} || X_i C12 + C12 Y_i ||_F^2

% return the fval and grad functions w.r.t. C12: a fMap from S1 to S2

Op1 = reshape(Op1, [], 1)';
Op2 = reshape(Op2, [], 1)';

if nargin < 3
    IfReversing = false;
end

if ~IfReversing             % do not revese the commutativity
    func_fval = @(C12) sum(cell2mat(...
        cellfun(@(X,Y) 0.5*norm(X*C12 - C12*Y,'fro')^2,...
        Op2, Op1, 'UniformOutput', false)...
        ), 2);
    func_grad = @(C12) sum(cell2mat(...
        cellfun(@(X,Y) reshape(X'*(X*C12 - C12*Y) - (X*C12 - C12*Y)*Y',[],1), ...
        Op2, Op1, 'UniformOutput', false)...
        ), 2);
else % used in finding symmetric maps
    func_fval = @(C12) sum(cell2mat(...
        cellfun(@(X,Y) 0.5*norm(X*C12 + C12*Y,'fro')^2,...
        Op2, Op1, 'UniformOutput', false)...
        ), 2);
    func_grad = @(C12) sum(cell2mat(...
        cellfun(@(X,Y) reshape(X'*(X*C12 + C12*Y) + (X*C12 + C12*Y)*Y',[],1), ...
        Op2, Op1, 'UniformOutput', false)...
        ), 2);
end
end