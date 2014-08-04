function Mshift = shift_matrix_Nd(mat, shift, dim)
% Mshift = shift_matrix_Nd(mat, shift, dim)
% INPUTS: 
% mat: matrix to which we're applying the shift
% shift: integer specifying the size of the shift 
% dim: dimension along which to shift the matrix
% OUTPUTS:
% Mshift: shifted matrix

if mod(shift,1) ~= 0
    error('Shift size must be integer valued');
end
mat_sz = size(mat);
mat_dim = length(mat_sz);
if mat_dim > 4
    error('Only support for up to 4d matrices');
end
if dim > mat_dim || dim < 1
    error('Invalid target dimension');
end
if shift > mat_sz(dim)
    error('Requested shift bigger than target dimension size');
end

%% for 0-shift, don't do anything
if shift == 0
    Mshift = mat;
    return;
end

%%
index_map = mod((0:(mat_sz(dim)-1))-shift,mat_sz(dim))+1;

if shift > 0
    zpad_inds = 1:shift;
elseif shift < 0
    zpad_inds = (mat_sz(dim) + shift + 1):mat_sz(dim);
end

switch mat_dim
    case 1
        Mshift = mat(index_map);
        Mshift(zpad_inds) = 0;
    case 2
        if dim == 1
            Mshift = mat(index_map,:);
            Mshift(zpad_inds,:) = 0;
        else
            Mshift = mat(:,index_map);
            Mshift(:,zpad_inds,:) = 0;
        end
    case 3
        if dim == 1
            Mshift = mat(index_map,:,:);
            Mshift(zpad_inds,:,:) = 0;
        elseif dim == 2
            Mshift = mat(:,index_map,:);
            Mshift(:,zpad_inds,:) = 0;
        else
            Mshift = mat(:,:,index_map);
            Mshift(:,:,zpad_inds) = 0;
        end
    case 4
        if dim == 1
            Mshift = mat(index_map,:,:,:);
            Mshift(zpad_inds,:,:,:) = 0;
        elseif dim == 2
            Mshift = mat(:,index_map,:,:);
            Mshift(:,zpad_inds,:,:) = 0;
        elseif dim == 3
            Mshift = mat(:,:,index_map,:);
            Mshift(:,:,zpad_inds,:) = 0;
        else
            Mshift = mat(:,:,:,index_map);
            Mshift(:,:,:,zpad_inds) = 0;
        end
end

