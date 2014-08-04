function [nks] = plot1dfilterbank_v2(ks,SDIM)
% USAGE: [nks] = plot1dfilterbank(ks,SDIM)
%   plot filterbank for 1D stimuli (i.e.) ala touryan or rust

flen = size(ks,1)/SDIM; 
nks  = size(ks,2); 
zrange = [min(ks(:)),max(ks(:))];  
% disp(zrange)
for ik = 1:nks; 
	subplot(nks,1,ik); plot2drf(ks(:,ik),SDIM,zrange);  
end

