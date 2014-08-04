function [beta] = predpow(psths,modpred)
%%USAGE: [beta] = predpow(psths,modpred)
%   Predictive Power a la Sahani LInden
 nreps    = size(psths,2); 
 meanpsth = mean(psths');  
 phat     = (nreps*var(meanpsth) - mean(var(psths)))/(nreps-1); 
 beta     = (var(meanpsth)-var(meanpsth-modpred))/phat; 
end

