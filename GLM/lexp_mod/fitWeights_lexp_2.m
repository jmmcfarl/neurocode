function fo = fitWeights_lexp_2(fo,g_mat,spkbs)
%[fo, LL] = fitWeights(fo,input,spkbs)
% fitWeights: fit modul-weights and intercept of the full GNLM

options.Display = 'iter';
options.maxFunEvals = 500;
options.maxIter = 200;
options.optTol = 1e-5;
options.progTol = 1e-7;

%compute binned spike vector
stimlen = size(g_mat,1);
Robs = zeros(1,stimlen);
ftable = tabulate(spkbs);
Robs(ftable(:,1)) = ftable(:,2);

nmods = length(fo.mods);

fX = [];
W0 = [];
min_x = min(fo.mods(1).nlx);
for i = 1:nmods
%     fX = [fX logexp(g_mat(:,i),fo.mods(i).beta)];
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
else
    error('Unsupported NL type');
end
    W0=[W0; fo.mods(i).w];
end
W0 = [W0; fo.const]; %tack on the constant term


lambda = ones(1,nmods)*fo.lambdaW;
lambda = [lambda 0]/sum(Robs);

[params LL] = L1General2_PSSas(@(K) lexp_weights_LLinternal(K, Robs, fX, fo),W0,lambda',options,0);


fo.const = params(end);
for imod = 1:nmods
    fo.mods(imod).w = params(imod);
end
fo.LL = LL;
