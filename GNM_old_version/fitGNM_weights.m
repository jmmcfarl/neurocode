function gnm = fitGNM_weights(gnm,g_mat,spkbs,silent,Pcon,l1_strength)

if nargin < 4
    silent = 1;
end
if nargin < 5
    Pcon = [];
end
if nargin < 6
    l1_strength = nan;
end
nmods = length(gnm.mods);

fX = [];
W0 = [];
min_x = min(gnm.mods(1).nlx);
for i = 1:nmods
    if strcmp(gnm.mods(i).nltype,'lexp')
    fX = [fX 1/gnm.mods(i).beta*log(1+exp(gnm.mods(i).beta*g_mat(:,i))) - 1/gnm.mods(i).beta*log(1+exp(gnm.mods(i).beta*min_x))];
    elseif strcmp(gnm.mods(i).nltype,'quad')
        fX = [fX g_mat(:,i).^2];
    elseif strcmp(gnm.mods(i).nltype,'lin')
        fX = [fX g_mat(:,i)];
    elseif strcmp(gnm.mods(i).nltype,'threshlin')
        temp = g_mat(:,i);
        temp(temp < 0) = 0;
        fX = [fX temp];
    elseif strcmp(gnm.mods(i).nltype,'rquad')
        temp = g_mat(:,i).^2;
        temp(g_mat(:,i) < 0) = 0;
        fX = [fX temp];
    elseif strcmp(gnm.mods(i).nltype,'uncon')
         fX = [fX nlin_proc_stim(g_mat(:,i),gnm.mods(i).nly,gnm.mods(i).nlx)];    
    else
        error('unsupported NL')
    end
    W0=[W0; gnm.mods(i).w];
end
W0 = [W0; gnm.spk_theta]; %tack on the constant term

%initialize parameters
lamrange = [];
lamrange2 = [];
Pmono = [];
hold_const = [];

%for L1 regularization of weights
if isnan(l1_strength)
    llist = [];
else
    llist = [l1_strength 1:(length(W0)-1)];
end
[fitp,grad] = GLMsolve_gnm(gnm, fX, spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const );
gnm.spk_theta = fitp.k(end);
% gnm.spk_theta = gnm.spk_theta - fitp.k(end);
for imod = 1:nmods
    gnm.mods(imod).w = fitp.k(imod);
end

%compute distribution of generating function
g = fX*fitp.k(1:end-1)+fitp.k(end);
[f,xi] = ksdensity(g);
gnm.fox = xi;
gnm.foy = f;
% fo.LL = fitp.LL;
% fo.LLnadj = fitp.LLnadj;

for imod =1:length(gnm.mods);
    if max(abs(gnm.mods(imod).k)) > 0
        cur_inputs = g_mat(:,imod);
        [foy,fox] = ksdensity(cur_inputs);
        gnm.mods(imod).fox = fox;
        gnm.mods(imod).foy = foy;
        
        [foys,foxs] = ksdensity(cur_inputs(spkbs));
        gnm.mods(imod).foxs = foxs;
        gnm.mods(imod).foys = foys;
    else
        gnm.mods(imod).fox = linspace(-3.1,3.1,10);
        gnm.mods(imod).foy = zeros(1,10);
        gnm.mods(imod).foxs = gnm.mods(imod).fox;
        gnm.mods(imod).foys = gnm.mods(imod).foy;
    end
end

end
