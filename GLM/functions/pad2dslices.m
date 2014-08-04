function [kcell] = pad2dslices(longk,spixids,sdim)
% USAGE: [kcell] = pad2dslices(longk,spixids,sdim)
%  takes the long k vector and converts it into a cell of matrices
nepix   = length(spixids); 
klen    = length(longk); 
nslices = klen/nepix; 

x       = flipud(reshape(longk,nslices,nepix)); 
kcell   = cell(nslices,1); 
for islice = 1:nslices; 
	kcell{islice} = padfiltermat(x(islice,:),spixids,sdim); 
end

