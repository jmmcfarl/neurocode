function nim = NIMadd_NLinput(nim,NLtype,NLsign,init_filt)
%
% nim = NIMadd_NLinput(nim,NLtype,NLsign,<init_filt>)
%
% Adds a subunit to the existing model nim.
%
% INPUTS:
%   nltype: specifies the type of upstream nonlinearity (can be 'lin','threshlin','quad', or 'nonpar')
%   nlsign: +1 or -1 to determine whether it's excitatory or suppressive
%   [init_filt]: vector specifying intial filter (default uses a random vector)
%
% OUTPUTS:
%   nim: New model struct
        
%%
kLen = length(nim.mods(1).filtK);
if nargin < 4
    init_filt = randn(kLen,1)/kLen * 0.1;
end

new_mod = nim.mods(1);
new_mod.NLtype = NLtype;
new_mod.filtK = init_filt;
new_mod.sign = NLsign;

nim.mods = cat(1,nim.mods,new_mod);