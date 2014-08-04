function [myglm] = fitGLM_temponly(myglm,T,spkbs,disp_type,its_per_round)

%INPUTS:
% myglm: model structure
% X: X-matrix
% spkbs: vector of spike bin indices
% disp_type: display type: 'tots','full','min','none'
% its_per_round: number of iterations to fit model before checking for user input to continue

optTol = 1e-4;
progTol = 1e-6;

nmods = length(myglm.mods);

if nargin < 6
    its_per_round = 200;
    check_user_in = 0;
else
    check_user_in = 1;
end

stimlen = size(T,1);

%initialize likelihood vecs
myglm.LP_seq = [];
myglm.LP_ax = []; %iteration number
myglm.LL_seq = [];

fprintf('INITIAL MODEL\n');
[ll0, ll0p] = getLLGLM_temponly(myglm,T,spkbs,disp_type);
myglm.LP = ll0p; myglm.LL = ll0;
myglm.LP_seq = [myglm.LP_seq ll0p];
myglm.LP_ax = [1];
myglm.LL_seq = [myglm.LL_seq ll0];

%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spkbs);
Robs(ftable(:,1)) = ftable(:,2);

eflag = 0;
n_opt_rounds = 0;
can_exit = 0;
while eflag == 0
    n_opt_rounds = n_opt_rounds + 1;
    fprintf('Optimization round %d\n',n_opt_rounds);
    [myglm eflag] = fitFULLRF_temponly(myglm,T,Robs,its_per_round,'lbgs',0,optTol,progTol);
    if check_user_in == 1
        R = input('Continue? (y/n)','s');
        if R == 'n'
            break
        end
    else
        break
    end
    eflag = 0;
end

%%
%if any subunits have all coefs to 0, move the 0 to the subunit weight
for imod = 1:nmods
    if all(myglm.mods(imod).k == 0)
        myglm.mods(imod).w = 0;
    end
end

[ll0, ll0p] = getLLGLM_temponly(myglm,T,spkbs,disp_type);
myglm.LP_seq(end) = ll0p;
myglm.LL_seq = [myglm.LL_seq ll0];

myglm.LL = ll0; myglm.LP = ll0p;


