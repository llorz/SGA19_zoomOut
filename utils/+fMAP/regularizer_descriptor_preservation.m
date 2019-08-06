function [func_fval, func_grad] = regularizer_descriptor_preservation(Fct1,Fct2)
% Preserve the projected descriptors:
%      @(C12) || C12 * Fct1 - Fct2   ||_F^2
% Fct_i: the projected (and normalized) corresponding descriptors of S_i
% return the fval and grad functions w.r.t. C12: a fMap from S1 to S2
func_fval = @(C12) 0.5*norm(C12*Fct1- Fct2,'fro')^2;
if nargout > 1
    func_grad = @(C12) reshape((C12*Fct1 - Fct2)*Fct1',[],1);
end
end