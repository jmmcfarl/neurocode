function fo = fitWeights_jmm(fo,fX,spkbs)
%[fo, LL] = fitWeights(fo,input,spkbs)
% fitWeights: fit modul-weights and intercept of the full GNLM

nmods = length(fo.mods); 

%compile all module weights into a single vector
W0 = [];
for i=1:nmods
    W0=[W0; fo.mods(i).w];
end
W0 = [W0; fo.const]; %tack on the constant term
 

%initialize parameters
silent = 0;
lamrange = [];
lamrange2 = [];
Pcon = [];
Pmono = [];
hold_const = [];
NLtype = 0;

%for L1 regularization of weights
llist = [fo.lambdaW];
llist = [llist 1:nmods]; 

[fitp,grad] = GLMsolve_jmm( fX, spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype );

fo.const = fitp.k(end);
for imod = 1:nmods
    fo.mods(imod).w = fitp.k(imod);
end
fo.LL = fitp.LL; 
fo.LLnadj = fitp.LLnadj;
end
