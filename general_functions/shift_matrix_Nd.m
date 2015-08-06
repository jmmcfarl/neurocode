function Mshift = shift_matrix_Nd(mat, shift, dim, pad_type)
% Mshift = shift_matrix_Nd(mat, shift, dim, pad_type)
% shift an Nd matrix 'mat' by 'shift' along dimension 'dim'. option for
% using zero of nan padding (default zero). Only works for up to 4d arrays
% at this point
% INPUTS:
%     mat: matrix to be shifted
%     shift: integer amount of shift
%     dim: dimension along which to shift
%     <pad_type>: 'zero' or 'nan', default 'zero', defines what kinds of padding to do to array
% OUTPUTS:
%     Mshift: shifted matrix

if nargin < 4
    pad_type = 'zero';
end

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
if strcmp(pad_type,'zero')
    pad_val = 0;
elseif strcmp(pad_type,'nan')
    pad_val = nan;
else
    error('Invalid padding type specified');
end

switch mat_dim
    case 1
        Mshift = mat(index_map);
        Mshift(zpad_inds) = pad_val;
    case 2
        if dim == 1
            Mshift = mat(index_map,:);
            Mshift(zpad_inds,:) = pad_val;
        else
            Mshift = mat(:,index_map);
            Mshift(:,zpad_inds) = pad_val;
        end
    case 3
        if dim == 1
            Mshift = mat(index_map,:,:);
                Mshift(zpad_inds,:,:) = pad_val;
        elseif dim == 2
            Mshift = mat(:,index_map,:);
            Mshift(:,zpad_inds,:) = pad_val;
        else
            Mshift = mat(:,:,index_map);
            Mshift(:,:,zpad_inds) = pad_val;
        end
    case 4
        if dim == 1
            Mshift = mat(index_map,:,:,:);
            Mshift(zpad_inds,:,:,:) = pad_val;
        elseif dim == 2
            Mshift = mat(:,index_map,:,:);
            Mshift(:,zpad_inds,:,:) = pad_val;
        elseif dim == 3
            Mshift = mat(:,:,index_map,:);
            Mshift(:,:,zpad_inds,:) = pad_val;
        else
            Mshift = mat(:,:,:,index_map);
            Mshift(:,:,:,zpad_inds) = pad_val;
        end
end

