function [pseudo_R2,mod_dev,null_dev] = pseudo_r2(r_obs,mod_r,nullmod_r)

if nargin < 3 || isempty(nullmod_r)
    nullmod_r = ones(size(r_obs))*mean(r_obs);
end

mod_LL = nansum(r_obs.*log(mod_r) - mod_r);
full_LL = nansum(r_obs.*log(r_obs) - r_obs);
mod_dev = 2*(full_LL-mod_LL);

null_mod_LL = nansum(r_obs.*log(nullmod_r) - nullmod_r);
null_dev = 2*(full_LL-null_mod_LL);

pseudo_R2 = 1-mod_dev./null_dev;


