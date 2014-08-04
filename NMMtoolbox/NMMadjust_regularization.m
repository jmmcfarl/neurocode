function nim = NMMadjust_regularization(nim,target_mods,varargin)
%
% nim = NMMadjust_regularization(nim,target_mods,'lambda_type1',lambda1_val,'lambda_type2',lambda2_val',...)
%
% Example: 
%   nim = NIMadjust_regularization(nim,[1 3 4],'lambda_d2T',100,'lambda_d2X',50);
%
% Takes a struct containing new regularization parameters to set for the
% vector of target subunits in the model nim
% 
% INPUTS:
%   nim: model struct
%   target_mods: vector of subunit indices to adjust regularization. Defaults to all mods.
%   Lambda_type: String that specifies type of regularization
%   Lambda_val: Value of regularization hyperparameters
% OUTPUTS:
%   nim: output model structure

%%
nmods = length(nim.mods);
if isempty(target_mods)
	target_mods = 1:nmods;
end

if any(~ismember(target_mods,1:nmods))
	error('Invalid target subunit');
end
if mod(length(varargin),2) ~= 0
	error('Invalid input format');
end

n_inputs = length(varargin)/2;
using_custom = false;
for ii = 1:n_inputs
	input_name{ii} = varargin{(ii-1)*2+1};
	if length(varargin{(ii-1)*2 + 2}) == 1
		input_val(ii,:) = varargin{(ii-1)*2 + 2} * ones(1,length(target_mods));
	else
		input_val(ii,:) = varargin{(ii-1)*2 + 2};
	end
		
	if ~isfield(nim.mods(1).reg_params,input_name{ii})
		error('Invalid regularization type');
	end
	if strcmp(input_name{ii},'lambda_custom')
		using_custom = true;
	end
end

for jj = 1:length(target_mods)
	tar = target_mods(jj);
	for ii = 1:n_inputs
		nim.mods(tar).reg_params = setfield(nim.mods(tar).reg_params,input_name{ii},input_val(ii,jj));
	end
end

if using_custom
    %fprintf('Using custom reg, setting other lambdas to 0 for target mods\n');
    for imod = target_mods
        nim.mods(imod).reg_params.lambda_dX = 0;
        nim.mods(imod).reg_params.lambda_dT = 0;
        nim.mods(imod).reg_params.lambda_d2XT = 0;
        nim.mods(imod).reg_params.lambda_d2X = 0;
        nim.mods(imod).reg_params.lambda_d2T = 0;
    end
end



