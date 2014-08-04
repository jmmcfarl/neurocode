function nim = NMMCinitialize_PSC( nim, target_mods, lambda_smooth, n_bfs, tent_spacing )
%
% Usage: nim = NMMCinitialize_PSC( nim, <target_mods>, <lambda_smooth>, <n_bfs>, <tent_spacing> )
%
% Initializes the specified model subunits to have nonparametric
% (tent-basis) upstream NLs
%
% INPUTS:
%     nim: input model structure
%     <target_mods>: vector specifying which subunits will be made to have nonparametric upstream NLs
%     <lambda_smooth>: specifies strength of smoothness regularization for the tent-basis coefs
%     <n_bfs>: Number of tent-basis functions to use
%     <tent_spacing>: 
% 
% OUTPUTS:
%     nim: output model

%%
Nmods = length(nim.mods);

if nargin < 2
	target_mods = 1:Nmods;
end
if nargin < 3
	lambda_smooth = 0;
end
if nargin < 4
	n_bfs = 10; % default number of tent bases
end
if nargin < 5
	tent_spacing = 1;
end

target_mods(target_mods <= 0) = [];

for imod = target_mods 	%nimr each module
    
	nim.mods(imod).PSC = zeros(n_bfs,1);
	nim.mods(imod).PSC(1) = 1;

	nim.mods(imod).PSCdt = tent_spacing;

end
