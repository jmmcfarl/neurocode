function nim = NMMCinitialize_Xconv( nim, target_mods, shift_range, lambda_smooth, lambda_L1 )
%
% Usage: nim = NMMCinitialize_Xconv( nim, <target_mods>, <shift_range>, <lambda_smooth>, <L1> )
%
% Initializes the specified model subunits to have nonparametric
% (tent-basis) upstream NLs
%
% Currently only works for 1-D convolutions
%
% INPUTS:
%     nim: input model structure
%     <target_mods>: vector specifying which subunits will be made to have nonparametric upstream NLs
%			<shift_range>: is +/- number given, but to Npix-1 in either side
%     <lambda_smooth>: specifies strength of smoothness regularization for the tent-basis coefs
%     <n_bfs>: Number of tent-basis functions to use
%     <tent_spacing>: 
% 
% OUTPUTS:
%     nim: output model

%%
Nmods = length(nim.mods);

if (nargin < 2) || isempty(target_mods)
	target_mods = 1:Nmods;
end
if (nargin < 3) || isempty(target_mods)
	shift_range = nim.stim_params.stim_dims(2)-1;
end
if (nargin < 4) || isempty(lambda_smooth)
	lambda_smooth = 0;
end
if nargin < 5
	lambda_L1 = 0;
end

target_mods(target_mods <= 0) = [];

for imod = target_mods 	%nimr each module
  
	nim.mods(imod).Xshifts = (-shift_range):shift_range;
	nim.mods(imod).Xconv = zeros(2*shift_range+1,1);
	nim.mods(imod).Xconv(shift_range+1) = 1;
	nim.mods(imod).reg_params.lambda_conv_L1 = lambda_L1; %second derivative of tent basis coefs on Tconv
    nim.mods(imod).reg_params.lambda_conv_d2X = lambda_smooth; %second derivative on Xconv

end

non_target_mods = setdiff(1:Nmods,target_mods);
for imod = non_target_mods
	nim.mods(imod).Xshifts = 0;
	nim.mods(imod).Xconv = 1;
	nim.mods(imod).reg_params.lambda_conv_L1 = 0; %second derivative of tent basis coefs on Tconv
    nim.mods(imod).reg_params.lambda_conv_d2X = 0; %second derivative on Xconv    
end
