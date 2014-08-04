function Xshifted = shift_mat_zpad(X,shift,dim)
%
% Xshifted = shift_mat_zpad(X,shift,dim)
%
% Takes a vector or matrix and shifts it along dimension dim by amount
% shift using zero-padding. Positive shifts move the matrix right or down

%%
sz = size(X);
if dim == 1
    if shift >= 0
        Xshifted = [zeros(shift,sz(2)); X(1:end-shift,:)];
    else
        Xshifted = [X(-shift+1:end,:); zeros(-shift,sz(2))];
    end
elseif dim == 2
    if shift >= 0
        Xshifted = [zeros(sz(1),shift) X(:,1:end-shift)];
    else
        Xshifted = [X(:,-shift:end); zeros(sz(1),-shift)];
    end
end

