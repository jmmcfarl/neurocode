function newK = filter_shift( k, dims, shift_amount, shift_dim )
%
% Usage: newK = filter_shift( k, dims, shift_amount, <shift_dim> )
%
% shift_dim defaults to first dim (max dimension is 3)

if nargin < 4
	shift_dim = 1;
end
if isempty(dims)
	% Assume one-dimensional 
	dims = size(k);
end

k3d = reshape( k, dims );
if shift_dim == 1
	if shift_amount > 0
		k3d((shift_amount+1):end,:,:) = k3d(1:end-shift_amount,:,:);
	else
		k3d(1:(end+shift_amount),:,:) = k3d((1-shift_amount):end,:,:);
	end
elseif shift_dim == 2
	if shift_amount > 0
		k3d(:,(shift_amount+1):end,:) = k3d(:,1:end-shift_amount,:);
	else
		k3d(:,1:(end+shift_amount),:) = k3d(:,(-shift_amount+1):end,:);
	end
else
	if shift_amount > 0
		k3d(:,:,(shift_amount+1):end) = k3d(:,:,1:end-shift_amount);
	else
		k3d(:,:,1:(end+shift_amount)) = k3d(:,:,(-shift_amount+1):end);
	end
end	

newK = reshape(k3d, prod(dims), 1, 1 );
