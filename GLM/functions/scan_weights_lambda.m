function [weight_mat,plls] = scan_weights_lambda(fo,stim,spkbs,lamvals)

nmods = length(fo.mods);

%initialize parameters
silent = 1;
lamrange = [];
lamrange2 = [];
Pcon = [];
Pmono = [];
hold_const = [];
NLtype = 0;

%compile a matrix representing the output of each NL module
fX    = [];
for i = 1:nmods
    fX = [fX procModul(fo.mods(i),stim)];
end

%compile all module weights into a single vector
W0 = [];
for i=1:nmods
    W0=[W0; fo.mods(i).w];
end
W0 = [W0; fo.const]; %tack on the constant term

    
%for L1 regularization of weights
llist = [fo.lambdaW];
llist = [llist 1:nmods];

%initialize weight_mat
weight_mat = nan(length(lamvals),nmods);
plls = nan(length(lamvals),1);

for i = 1:length(lamvals)
    
    fprintf('Trying %d of %d lambda values\n',i,length(lamvals));
    llist(1) = lamvals(i);
    [fitp,grad] = GLMsolve_jmm( fX, spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype );
    weight_mat(i,:) = fitp.k(1:end-1);
    plls(i) = fitp.LL;
    
end