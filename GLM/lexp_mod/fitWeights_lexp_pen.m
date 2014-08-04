function fo = fitWeights_lexp_pen(fo,X,spkbs,silent,Pcon)
%[fo, LL] = fitWeights(fo,input,spkbs)
% fitWeights: fit modul-weights and intercept of the full GNLM

if nargin < 4
    silent = 1;
end
if nargin < 5
    Pcon = [];
end
nmods = length(fo.mods);

k_mat = get_k_mat(fo);
g_mat = X*k_mat;

fX = [];
W0 = [];
min_x = min(fo.mods(1).nlx);
for i = 1:nmods
    if strcmp(fo.mods(i).nltype,'lexp')
    fX = [fX 1/fo.mods(i).beta*log(1+exp(fo.mods(i).beta*g_mat(:,i))) - 1/fo.mods(i).beta*log(1+exp(fo.mods(i).beta*min_x))];
    elseif strcmp(fo.mods(i).nltype,'quad')
        fX = [fX g_mat(:,i).^2];
    elseif strcmp(fo.mods(i).nltype,'lin')
        fX = [fX g_mat(:,i)];
    elseif strcmp(fo.mods(i).nltype,'threshlin')
        temp = g_mat(:,i);
        temp(temp < 0) = 0;
        fX = [fX temp];
    elseif strcmp(fo.mods(i).nltype,'rquad')
        temp = g_mat(:,i).^2;
        temp(g_mat(:,i) < 0) = 0;
        fX = [fX temp];
    elseif strcmp(fo.mods(i).nltype,'uncon')
         fX = [fX nlin_proc_stim(g_mat(:,i),fo.mods(i).nly,fo.mods(i).nlx)];    
    else
        error('unsupported NL')
    end
    W0=[W0; fo.mods(i).w];
end
W0 = [W0; fo.const]; %tack on the constant term

%initialize parameters
lamrange = [];
lamrange2 = [];
Pmono = [];
hold_const = [];
if strcmp(fo.spk_nl,'logexp')
    NLtype = 0;
elseif strcmp(fo.spk_nl,'exp')
    NLtype = 1;
end

[nll, pnll,lpen] = getLLGLM_lexp(fo,X,spkbs,'none');
tot_pens = lpen.l2x + lpen.ksmooth + lpen.lapl + lpen.laplXT + lpen.loc + lpen.l1x';
lambda_Ws = zeros(nmods,1);
for i = 1:nmods
    if strcmp(fo.mods(i).nltype,'lin') || strcmp(fo.mods(i).nltype,'threshlin')
        lambda_Ws(i) = tot_pens(i)
    elseif strcmp(fo.mods(i).nltype,'quad') || strcmp(fo.mods(i).nltype,'rquad')
        lambda_Ws(i) = tot_pens(i)
    else
        error('Cant penalize weights with this internal NL')
    end
end

%for L1 regularization of weights
llist = [];
if fo.lambdaW > 0
    llist = [fo.lambdaW];
    llist = [llist 1:nmods]; %all constrained to be positive
end

[fitp,grad] = GLMsolve_jmm( fX, spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype );

fo.const = fitp.k(end);
for imod = 1:nmods
    fo.mods(imod).w = fitp.k(imod);
end

%compute distribution of generating function
g = fX*fitp.k(1:end-1)+fitp.k(end);
[f,xi] = ksdensity(g);
fo.fox = xi;
fo.foy = f;
% fo.LL = fitp.LL;
% fo.LLnadj = fitp.LLnadj;

for imod =1:length(fo.mods);
    if max(abs(fo.mods(imod).k)) > 0
        cur_inputs = g_mat(:,imod);
        [foy,fox] = ksdensity(cur_inputs);
        fo.mods(imod).fox = fox;
        fo.mods(imod).foy = foy;
        
        [foys,foxs] = ksdensity(cur_inputs(spkbs));
        fo.mods(imod).foxs = foxs;
        fo.mods(imod).foys = foys;
    else
        fo.mods(imod).fox = linspace(-3.1,3.1,10);
        fo.mods(imod).foy = zeros(1,10);
        fo.mods(imod).foxs = fo.mods(imod).fox;
        fo.mods(imod).foys = fo.mods(imod).foy;
    end
end

end
