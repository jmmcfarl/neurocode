function new_X = tb_proc_stim(X,tent_space,flen)

if nargin < 3
    flen = 1;
end

%create a tent-basis (triangle) filter
tent_filter = [(1:tent_space)/tent_space 1-(1:tent_space-1)/tent_space]/tent_space;

nPix = size(X,2)/flen;
X = reshape(X,size(X,1),flen,nPix);

%apply to the stimulus
new_X = zeros(size(X));
for i = 1:length(tent_filter)
    new_X = new_X + shift_matrix_Nd(X,i-tent_space,3)*tent_filter(i);
end
new_X = new_X(:,:,1:tent_space:end);
new_X = reshape(new_X,size(new_X,1),[]);
