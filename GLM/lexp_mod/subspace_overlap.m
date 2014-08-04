function O = subspace_overlap(A,B)

%computes overlap between subspaces defined by column vectors of A an B.
%dim(A) must be >= dim(B)

if size(A,2) == size(B,2)
    k = size(A,2);
    Aorth = gramschmidtorth(A);
    Borth = gramschmidtorth(B);
    O = abs(det(Aorth'*Borth))^(1/k)/(abs(det(Aorth'*Aorth))^(1/(2*k)) ...
        * abs(det(Borth'*Borth))^(1/(2*k)));
elseif size(A,2) > size(B,2)
    k = size(B,2);
    Aprojmat = A/(A'*A)*A';
    B_A = Aprojmat*B;
    B_A_orth = gramschmidtorth(B_A);
    B_orth = gramschmidtorth(B);
    O = abs(det(B_A_orth'*B_orth))^(1/k)/(abs(det(B_A_orth'*B_A_orth))^(1/(2*k)) ...
        * abs(det(B_orth'*B_orth))^(1/(2*k)));    
else
    error('dim(A) must be >= dim(B)');
end


