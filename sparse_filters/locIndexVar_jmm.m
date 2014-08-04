function [locindex] = locIndexVar_jmm(tk,rfcenter,SDIM)
% USAGE:  [locindex] = locIndexVar(basevecs,rfcenter,beta,SDIM)
%   localization index reminding of computation of variance -- squared deviation
%   from rfcenter and rfshape is interpreted similar to a probability distribution
%   but not normalized (i.e. sum(vars) does not equal 1)

flen     = length(tk)/SDIM; 
vars     = var(reshape(l2norm(tk),flen,SDIM)); %temporal variance in each dimension of normalized filter

distvec  = ((1:SDIM)-rfcenter); %distance from the RF center point
locindex = (distvec.*distvec)*vars'; %localization measure is the temporal variance at each pixel weighted by the squared distance from the RF center 
%locindex = 10*abs(distvec)*vars';  
end

