function [func_fval, func_grad] = regularizer_laplacian_commutativity(Ev1, Ev2)
% standard Laplacian commutativity term:
%               @(C12) || C12 diag(Ev1) - diag(Ev2) C12 ||_F^2
% Ev_i: the eigenvalues of shape S_i
% return the fval and grad functions w.r.t. C12: a fMap from S1 to S2
k1 = length(Ev1);
k2 = length(Ev2);

Ev1 = reshape(Ev1, [], 1);
Ev2 = reshape(Ev2, [], 1);

% Dlb = (repmat(Ev2/sum(Ev2)*sum(Ev1), [1,k1]) - repmat(Ev1'/sum(Ev1)*sum(Ev2), [k2,1])).^2;
Dlb = (repmat(Ev2, [1,k1]) - ...
    repmat(Ev1', [k2,1])).^2;
Dlb = Dlb/norm(Dlb, 'fro');

func_fval = @(C12) sum(sum((C12.^2 .* Dlb)/2));
if nargout > 1
    func_grad = @(C12) reshape(C12.*Dlb,[],1);
end

end