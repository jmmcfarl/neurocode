function nim = NMMadd_NLinput( nim, NLtype, NLsign, con_type, Xtarget, init_filt )
%
% Usage: nim = NMMadd_NLinput( nim, NLtype, NLsign, <con_type>, <Xtarget>, <init_filt> )
%
% Adds a subunit to the existing model nim.
%
% INPUTS:
%   nltype: specifies the type of upstream nonlinearity (can be 'lin','threshlin','quad', or 'nonpar')
%   nlsign: +1 or -1 to determine whether it's excitatory or suppressive
%   <con_type>: 0 (default) always positive, monotonic increaseing, first element anchor at 0
%               1 all between 0 and 1, monotonic decreasing, no anchor
%               2 all between 0 and 1, monotonic decreasing, first element anchor at 1
%               x (>2) all between 0 and 1, no constraint, anchor element set here
%   Xtarget: list of which X-matrices to operate on
%   stim_params: leave blank if Xtarget is same, otherwise will add new Xta
%   [init_filt]: vector specifying intial filter (default uses a random vector)
%
% OUTPUTS:
%   nim: New model struct
        
%%
kLen = length(nim.mods(1).filtK);
if nargin < 6
	init_filt = randn(kLen,1)/kLen * 0.1;
end
if (nargin < 5) || isempty(Xtarget)
	Xtarget = 1;
end
if (nargin < 4) || isempty(con_type)
	con_type = 0;
end

if Xtarget > length(nim.stim_params)
	error( 'Need to define more stim parameters before adding new X-target.' )
end

new_mod = nim.mods(1);
new_mod.Xtarget = Xtarget;
new_mod.NLtype = NLtype;
new_mod.con_type = con_type;
new_mod.filtK = init_filt;
new_mod.sign = NLsign;
new_mod.Kcon = 0;

nim.mods = cat(1,nim.mods,new_mod);