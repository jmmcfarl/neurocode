function [bestk,LL] = fitGLM_weights_gabor(X,spkbs,init_k,silent)

[NT,k] = size(X);

% create lamrange matrix
lnl2 = 500;
lamrange2 = [lnl2 1 k];

% nlcon = 1;
% Pcon = [nlcon*(1:k)];
Pcon = [];
llist = [];
Pmono = [];
Pzero = [];

init_k = [init_k; 0]; %tack on the constant term

hold_const = [];

%%
Ncon = length(Pcon);

if Ncon > 0
    mults2 = ones(Ncon,1);
    mults2(Pcon < 0) = -1;
    Pcon = abs(Pcon);
end

L = length(init_k);
A = zeros(Ncon,L);
B = zeros(Ncon,1);if silent == 0
    opts = optimset('GradObj','on','LargeScale','off','Display','iter','MaxIter',2000,'MaxFunEvals',10000,'TolFun',1e-7);
    %     opts = optimset('GradObj','on','Algorithm','active-set','Display','iter','MaxIter',200,'MaxFunEvals',10000,'TolFun',TolFun);
else
    opts = optimset('GradObj','on','LargeScale','off','MaxIter',2000,'Display','off','MaxFunEvals',10000,'TolFun',1e-7);
    %     opts = optimset('GradObj','on','Algorithm','active-set','MaxIter',200,'Display','off','MaxFunEvals',10000,'TolFun',TolFun);
end


for i = 1:Ncon
    A(i,Pcon(i)) = -mults2(i);
end

%%
if silent == 0
    opts = optimset('GradObj','on','LargeScale','off','Display','iter','MaxIter',2000,'MaxFunEvals',10000,'TolFun',1e-7);
    %     opts = optimset('GradObj','on','Algorithm','active-set','Display','iter','MaxIter',200,'MaxFunEvals',10000,'TolFun',TolFun);
else
    opts = optimset('GradObj','on','LargeScale','off','MaxIter',2000,'Display','off','MaxFunEvals',10000,'TolFun',1e-7);
    %     opts = optimset('GradObj','on','Algorithm','active-set','MaxIter',200,'Display','off','MaxFunEvals',10000,'TolFun',TolFun);
end

NLtype = 1;
S = X(spkbs,:);
if Ncon > 0
[bestk,LL,eflag,outp,lam,grad] = fmincon( @(K) LLexp_jmm(K,X,S,[],lamrange2,llist,hold_const), init_k, A, B, [], [],[],[],[], opts );
else
 [bestk,LL,eflag,outp,lam,grad] = fminunc( @(K) LLexp_jmm(K,X,S,[],lamrange2,llist,hold_const), init_k, opts );   
end
%%
