function [nll, pnll,lpen] = getLL_jmm_STCBF(kx,glmod,spkbs)
%% computes likelihood of generator potential
%% uses penalty terms for both nl and history parameters 

%% pure ll based on model prediction (e.g. likelihood)
kx(kx > 50)    = 50; %saturate input to spiking NL
if strcmp(glmod.spk_nl,'logexp')
    rx             = log(1+exp(kx)); %apply spiking NL tp get predicted rate
elseif strcmp(glmod.spk_nl,'exp')
    rx = exp(kx);
end
rs             = rx(spkbs); %predicted rate at spike times                  
rs(rs < 1e-10) = 1e-10; %impose minimum on predicted rate

%% likelihood given by point process model
ll = sum(log(rs)) - sum(rx);

%% regularisaton penalty (e.g. enforcing smth like posterior)
lambdaW = glmod.lambdaW;
sdim = glmod.mods(1).SDIM;
nmods = length(glmod.mods);

nspks = length(spkbs);

%initialize penalty terms
lpen.h = zeros(nmods,1);
lpen.h2 = zeros(nmods,1);
lpen.nl = zeros(nmods,1);
lpen.nl2 = zeros(nmods,1);
lpen.loc = zeros(nmods,1);
lpen.w = zeros(nmods,1);
for i = 1:length(glmod.mods)
  mod     = glmod.mods(i); %for current module
  nlx = mod.nlx;
  nly = mod.nly;
  
  lpen.h(i)  = (mod.lh  * ssq(diff(mod.h)))/nspks; %penalty on PSC slope = lambda_h * (total slope) 
%   lpen.nl(i) = (mod.lnl * ssq(diff(mod.nly)))/nspks; %penalty on NL slope = lambda_NL * (total slope) (not true slope though!!)
  lpen.nl(i) = (mod.lnl * ssq(diff(nly)./diff(nlx)))/nspks; %penalty on NL slope = lambda_NL * (total slope) (not true slope though!!)
  
  lpen.h2(i) = (mod.lh2 * sum((2*mod.h(2:end-1) - mod.h(1:end-2) - mod.h(3:end)).^2))/nspks;
%   lpen.nl2(i) = (mod.lnl2 * sum((2*mod.nly(2:end-1) - mod.nly(1:end-2) - mod.nly(3:end)).^2))/nspks; %this is not the true slope calc!)
 
dy_dx2 = (nly(3:end)-nly(2:end-1))./(nlx(3:end)-nlx(2:end-1));
dy_dx1 = (nly(2:end-1)-nly(1:end-2))./(nlx(2:end-1)-nlx(1:end-2));
nlx_mids = (nlx(2:end)+nlx(1:end-1))/2;
d2y_d2x = (dy_dx2-dy_dx1)./diff(nlx_mids);

% lpen.nl2(i) = (mod.lnl2 * sum((2*mod.nly(2:end-1) - mod.nly(1:end-2) - mod.nly(3:end)).^2))/nspks; %this is not the true slope calc!)
lpen.nl2(i) = (mod.lnl2 * sum(d2y_d2x.^2))/nspks; %this is not the true slope calc!)

lpen.w(i) = lambdaW*abs(mod.w)/nspks;

  %kernel localization penalty
  lpen.loc(i) = mod.locLambda*kernel_std(glmod.STCbasis,mod.STCcf,sdim)/nspks;
%     pll = pll - glmod.mods(i).locLambda*kernel_entropy(glmod.STCbasis,mod.STCcf,sdim);

%for prior covmat penalty on k's
%     cur_kern = glmod.STCbasis*glmod.mods(i).STCcf;
%     pll = pll - cur_kern'*glmod.mods(i).prior_precisionmat*cur_kern;

end


%% penalty terms for spiketrain moduls?
nll = -ll/nspks;
pnll = nll + sum(lpen.h) + sum(lpen.nl) + sum(lpen.h2) + sum(lpen.nl2) + sum(lpen.w) + sum(lpen.loc); 

end