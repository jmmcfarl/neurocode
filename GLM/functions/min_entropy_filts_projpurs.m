function [loc_kerns loc_cfs] = min_entropy_filts_projpurs(STCbvs,ncomps,sdim)
% USAGE: [subunits rflocs] = findLocalSubunits(nestcs,ncomps,SDIM)

[flen,nbvs] = size(STCbvs);

loc_kerns = zeros(flen,ncomps); %initialize localized subunits
loc_cfs = zeros(nbvs,ncomps);

ooptions = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-5,'TolX',1e-5);   

beta0    = ones(nbvs,1)/sqrt(nbvs); %uniform intial guess

cur_STC_space = STCbvs;

%loop over localized filters
for icomp = 1:ncomps; 
    
    %search for linear combo of STC filters which maximize the power
    %concentration measure
	[argmin,minval] = fminsearch(@(x) kernel_entropy(cur_STC_space,x,sdim),beta0,ooptions);
    
	betas  = l2norm(argmin); %normalize the linear combo
	trf    = l2norm(cur_STC_space*betas); %normalized projection onto first component 
    
    loc_kerns(:,icomp) = trf;
        
	cur_STC_space = cur_STC_space - trf*(cur_STC_space'*trf)'; %subtract off projection of STC matrix onto loc comp
	cur_STC_space = cur_STC_space*diag(1./sqrt(diag(cur_STC_space'*cur_STC_space))); %rescale to have equal variance in each remaining direction
    
	disp(sprintf('--- subunit %i: %f %f',icomp,minval));
end 

loc_cfs = (loc_kerns'*STCbvs)';