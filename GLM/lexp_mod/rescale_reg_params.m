function [new_mod,scale] = rescale_reg_params(f0,X,defmod)

new_mod = f0;
g_mat = X*get_k_mat(f0);
nmods = length(f0.mods);
scale = zeros(nmods,1);
for imod = 1:nmods
    mod = f0.mods(imod);
    if strcmp(mod.nltype,'lexp')
        fg = 1/mod.beta*log(1+exp(mod.beta*(g_mat(:,imod)-mod.theta))) - 1/mod.beta*log(1+exp(mod.beta*(min_x-mod.theta)));
    elseif strcmp(mod.nltype,'quad')
        fg = g_mat(:,imod).^2;
    elseif strcmp(mod.nltype,'lin')
        fg = g_mat(:,imod);
    elseif strcmp(mod.nltype,'threshlin')
        fg = g_mat(:,imod);
        fg(fg < 0) = 0;
    elseif strcmp(mod.nltype,'rquad')
        fg = g_mat(:,imod).^2;
        fg(g_mat(:,imod) < 0) = 0;
    elseif strcmp(mod.nltype,'uncon')
        fg    = nlin_proc_stim(g_mat(:,imod),mod.nly,mod.nlx); %pass internal output through NL
    else
        error('Invalid nl type');
    end
    fg = fg*mod.w;
    scale(imod) = std(fg);
    mod.lambda_dX = defmod.lambda_dX/scale;
    mod.lambda_dT = defmod.lambda_dT/scale;
    mod.lambda_L1x = defmod.lambda_L1x/scale;
    mod.lambda_d2X = defmod.lambda_d2X/scale;
    mod.lambda_d2XT = defmod.lambda_d2XT/scale;
%     mod.lambda_L2X = defmod.lambda_L2X/scale;
    mod.locLambda = defmod.locLambda/scale;
end