function penalties = get_NMM_reg_penalties( nim, regmat_custom)
%

%% Set parameters
if nargin < 2
    regmat_custom = [];
end
Nmods = length(nim.mods);
%% CREATE L2 REGULARIZATION MATRICES
L2_mats = create_L2_matrices_NMM(nim);
L2_mats.custom = regmat_custom;

%% COMPUTE L2 PENALTIES
smooth_penalty = zeros(Nmods,1);
deriv_penalty = zeros(Nmods,1);
ridge_penalty = zeros(Nmods,1);
sparse_penalty = zeros(Nmods,1);
custom_penalty = zeros(Nmods,1);
for n = 1:Nmods
    cur_kern = nim.mods(n).filtK;
    if nim.mods(n).reg_params.lambda_dT > 0
        deriv_penalty(n) = deriv_penalty(n) + nim.mods(n).reg_params.lambda_dT*sum((L2_mats.L2_dT{nim.mods(n).Xtarget} * cur_kern).^2);
    end
    if nim.mods(n).reg_params.lambda_d2T > 0
        smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2T*sum((L2_mats.L2_d2T{nim.mods(n).Xtarget} * cur_kern).^2);
    end
%     if nim.stim_params.stim_dims(2) > 1 && nim.stim_params.stim_dims(3) == 1
    if nim.stim_params(nim.mods(n).Xtarget).stim_dims(2) > 1 %if spatial dims
        if nim.mods(n).reg_params.lambda_dX > 0
            deriv_penalty(n) = deriv_penalty(n) + nim.mods(n).reg_params.lambda_dX*sum((L2_mats.L2_dX{nim.mods(n).Xtarget} * cur_kern).^2);
        end
        if nim.mods(n).reg_params.lambda_d2X > 0
            smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2X*sum((L2_mats.L2_d2X{nim.mods(n).Xtarget} * cur_kern).^2);
        end
        if nim.mods(n).reg_params.lambda_d2XT > 0
            smooth_penalty(n) = smooth_penalty(n) + nim.mods(n).reg_params.lambda_d2XT*sum((L2_mats.L2_d2XT{nim.mods(n).Xtarget} * cur_kern).^2);
        end
    end
    if nim.mods(n).reg_params.lambda_L2 > 0
        ridge_penalty(n) = nim.mods(n).reg_params.lambda_L2*(cur_kern' * cur_kern);
    end
    if nim.mods(n).reg_params.lambda_L1 > 0
        sparse_penalty(n) = nim.mods(n).reg_params.lambda_L1*sum(abs(cur_kern));
    end
    %for custom regularization
    if nim.mods(n).reg_params.lambda_custom > 0
			if isempty(L2_mats.custom)
% 				disp('Warning: penalty calculation is off because L2_mats.custom is not included.')
			else
        custom_penalty(n) = custom_penalty(n) + nim.mods(n).reg_params.lambda_custom*sum((L2_mats.custom * cur_kern).^2);
			end
    end
    
end

%penalized LL
total_penalty = sum(smooth_penalty) + sum(ridge_penalty) + sum(deriv_penalty) + sum(sparse_penalty) + sum(custom_penalty);

%%
penalties.smoothness = smooth_penalty;
penalties.deriv = deriv_penalty;
penalties.ridge = ridge_penalty;
penalties.sparse = sparse_penalty;
penalties.custom = custom_penalty;