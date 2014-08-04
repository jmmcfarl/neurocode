function [pseudo_R2,mod_dev,null_dev] = pseudo_r2(r_obs,mod_r,nullmod_r)
% [pseudo_R2,mod_dev,null_dev] = pseudo_r2(r_obs,mod_r,[nullmod_r])
% INPUTS:
% r_obs: observed spike count vector
% mod_r: model-predicted firing rate
% [nullmod_r]: firing rate prediction of the null model (assumed constant if none given).
% OUTPUTS:
% pseudo_R2: Estimated pseudo-R2
% mod_dev: Deviance of the model
% null_dev: Deviance of the null model

%if no prediction specified for the null model, use constant rate
%prediction
if nargin < 3 || isempty(nullmod_r)
    nullmod_r = ones(size(r_obs))*mean(r_obs);
end

mod_LL = nansum(r_obs.*log(mod_r) - mod_r);
full_LL = nansum(r_obs.*log(r_obs) - r_obs);
mod_dev = 2*(full_LL-mod_LL);

null_mod_LL = nansum(r_obs.*log(nullmod_r) - nullmod_r);
null_dev = 2*(full_LL-null_mod_LL);

pseudo_R2 = 1-mod_dev./null_dev;


