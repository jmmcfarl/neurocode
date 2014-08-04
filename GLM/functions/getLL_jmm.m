function [nll, pnll] = getLL_jmm(kx,glmod,spkbs)
%% computes likelihood of generator potential
%% uses penalty terms for both nl and history parameters 

%% pure ll based on model prediction (e.g. likelihood)
kx(kx > 20)    = 20; %saturate input to spiking NL
rx             = log(1+exp(kx)); %apply spiking NL tp get predicted rate
rs             = rx(spkbs); %predicted rate at spike times                  
rs(rs < 1e-10) = 1e-10; %impose minimum on predicted rate

%% likelihood given by point process model
ll = sum(log(rs)) - sum(rx);

%% regularisaton penalty (e.g. enforcing smth like posterior)
lambdaW = glmod.lambdaW;
pll = ll; %jmm 11/1/11 added
for i = 1:length(glmod.mods)
  mod     = glmod.mods(i); %for current module
  ssqs.h  = mod.lh  * ssq(diff(mod.h)); %penalty on PSC slope = lambda_h * (total slope) 
  ssqs.nl = mod.lnl * ssq(diff(mod.nly)); %penalty on NL slope = lambda_NL * (total slope)
 
  %shouldn't this be initialized and decremented each time to get the total penalized LL across modules?
  pll     = pll - ssqs.h - ssqs.nl; %jmm 11/1/11, changed from pll = ll - ssqs.h - ssqs.nl;
  pll     = pll - lambdaW*abs(mod.w);
end


%% penalty terms for spiketrain moduls?


nll  = - ll / length(spkbs); %raw neg LL per spike
pnll = - pll / length(spkbs); %penalized neg LL per spike
end