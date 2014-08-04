function nim_out = NMMset_regpenalties(nim,reg_type,pen_vals,targets)
%
%

%%
nim_out = NMMadjust_regularization(nim,targets,reg_type,ones(length(targets),1));
L2_mats = create_L2_matrices_NMM(nim_out);
if strcmp('reg_type','lambda_custom')
    L2_mats.custom = regmat_custom;
end
%%
Nmods = length(targets);
reg_norm = zeros(Nmods,1);
for n = 1:Nmods
    cur_kern = nim.mods(targets(n)).filtK;
    if strcmp(reg_type,'lambda_dT')
        reg_norm(n) = sum((L2_mats.L2_dT{nim.mods(targets(n)).Xtarget} * cur_kern).^2);
    end
    if strcmp(reg_type,'lambda_d2T')
        reg_norm(n) = sum((L2_mats.L2_d2T{nim.mods(targets(n)).Xtarget} * cur_kern).^2);
    end
    if strcmp(reg_type,'lambda_dX')
        reg_norm(n) = sum((L2_mats.L2_dX{nim.mods(targets(n)).Xtarget} * cur_kern).^2);
    end
    if strcmp(reg_type,'lambda_d2X')
        reg_norm(n) = sum((L2_mats.L2_d2X{nim.mods(targets(n)).Xtarget} * cur_kern).^2);
    end
    if strcmp(reg_type,'lambda_d2XT')
        reg_norm(n) = sum((L2_mats.L2_d2XT{nim.mods(targets(n)).Xtarget} * cur_kern).^2);
    end
    if strcmp(reg_type,'lambda_L2')
        reg_norm(n) = (cur_kern' * cur_kern);
    end
    if strcmp(reg_type,'lambda_L1')
        reg_norm(n) = sum(cur_kern);
    end
    %for custom regularization
    if strcmp(reg_type,'lambda_custom')
        reg_norm(n) = sum((L2_mats.custom * cur_kern).^2);
    end
end

%%
new_lambdas = pen_vals(:)./reg_norm(:);
nim_out = NMMadjust_regularization(nim_out,targets,reg_type,new_lambdas);
