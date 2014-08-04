function [Ov,vec_var] = subspace_overlap2(A,B)

B = bsxfun(@rdivide,B,sqrt(sum(B.^2)));

Aprojmat = A/(A'*A)*A';
B_A = Aprojmat*B;

vec_var = sum(B_A.^2);
Ov = mean(vec_var);

