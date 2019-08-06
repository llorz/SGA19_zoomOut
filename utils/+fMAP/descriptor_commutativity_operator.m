function CommOp = descriptor_commutativity_operator(desc, B, A)
% for each given descriptor, construct a corresponding operator to be
% preserved by the fMap via commutativity
% desc: a set of descriptors
% B: the basis where the fMap lives in
% A: the area matrix of the corresponding shape

num_desc = size(desc,2);
num_basis = size(B,2);
CommOp = cell(num_desc,1);
for i = 1:num_desc
    CommOp{i} = B'*A*(repmat(desc(:,i), [1,num_basis]).*B);
end

end