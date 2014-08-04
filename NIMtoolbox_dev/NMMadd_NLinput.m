function nim = NMMadd_NLinput( nim, NLtype, NLsign, Xtarget, init_filt, reg_params )
%
% Usage: nim = NMMadd_NLinput( nim, NLtype, NLsign, <Xtarget>, <init_filt>,<reg_params> )
%
% Adds a subunit to the existing model nim.
%
% INPUTS:
%   nltype: specifies the type of upstream nonlinearity (can be 'lin','threshlin','quad', or 'nonpar')
%   nlsign: +1 or -1 to determine whether it's excitatory or suppressive
%   Xtarget: list of which X-matrices to operate on
%   [init_filt]: vector specifying intial filter (default uses a random vector)
%   [reg_params]: optional specify reg params struct for the added input
%
% OUTPUTS:
%   nim: New model struct
        
%%
%%
if (nargin < 4) || isempty(Xtarget)
	Xtarget = 1;
end
if Xtarget > length(nim.stim_params)
	error( 'Need to define more stim parameters before adding new X-target.' )
end

kLen = prod(nim.stim_params(Xtarget).stim_dims);
if nargin < 5
	init_filt = randn(kLen,1)/kLen * 0.1;
end
if nargin < 6
    reg_params = [];
end
if Xtarget > length(nim.stim_params)
	error( 'Need to define more stim parameters before adding new X-target.' )
end

new_mod = nim.mods(1);
new_mod.Xtarget = Xtarget;
new_mod.NLtype = NLtype;
new_mod.filtK = init_filt;
new_mod.sign = NLsign;
if ~isempty(reg_params)
    new_mod.reg_params = reg_params;
end

nim.mods = nim.mods(:);
nim.mods = cat(1,nim.mods,new_mod);