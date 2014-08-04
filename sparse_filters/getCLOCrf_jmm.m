function [tlrf] = getCLOCrf_jmm(tpos,tks,SDIM,ninits,ooptions)
%USAGE: [tlrf] = getCLOCrf(tpos,tks,SDIM,ninits)
%   find most localized filter centered on position tpos using ninits

nbvs  = size(tks,2); %number of dimensions in the STC space
flen  = size(tks,1)/SDIM; 
lrfs  = zeros(SDIM*flen,ninits); %initialize localized filters
locis = zeros(ninits,1); 
for irep = 1:ninits; %for ninits initializations
	beta0 = l2norm(randn(nbvs,1)); %initialize guess to be normalized gaussian noise
	[targmin,minval] = fminsearch(@(x) locIndexVar_jmm(tks*x,tpos,SDIM),beta0,ooptions); %find linear combo of STCs producing the most localized filter
	lrfs(:,irep)     = tks*l2norm(targmin); %compute and store the resultant filter
	locis(irep)      = minval; %store the localization value
end;
disp(locis');
[lmin,lmindex] = min(locis); %take the most localized filter across the different initializations
tlrf = l2norm(lrfs(:,lmindex));
