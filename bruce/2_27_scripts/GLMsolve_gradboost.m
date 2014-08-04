function [fitp,grad] = GLMsolve_gradboost( X, cur_mod_out, spksN, K0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, TolFun, nlxrange )
%
% Usage: [fitp,grad] = GLMsolve( X, spksN, K0, <silent>, <lamrange>, <lamrange2>, <Pcon>, <Pmono>, <llist>, <hold_const>, <NLtype> )
%
% 'X' represents the data as a function of time used to predict the spike response 'spksN'
%    -- the output of the linear stage: g = X . k + b then a function of the vector of
%       parameters k and the offset b
% 'K0' is an initial guess of [k b]
% 'silent' > 0 suppresses the output of matlabs gradient descent
% 'lamrange' is a Nx3 matrix that specifies slope-penalty terms (ridge regression)
%     in form [multiplier range_beg range_end]
%     e.g. [0.1 1 12; 1 13 24;] applies penalty to LL of -0.1*(diff(k(1:12))^2 - 1*(diff(k(13:24)))^2
% 'Pcon' is a list of parameters of k that are constrained to be only positive or only negative
%     its comprisef of list of k-indices, which are positive if +constr, and negative if -constr
%     e.g. Pcon = [10 11 12 13 -14 -15 -16 -17] constrains k(10:13) > 0 and k(14:17) < 0
% 'Pmono' is a Mx2 list of sections of the space that should be monotonically increasing or decreasing
%     e.g. Pmon = [1 10; -11 -20;] constrains k(1:10) to be monot increasing and k(11:20) to be m.dec.
% 'llist' applies a sparseness penalty term of form lambda*abs(k).
%       llist is of the form [Lambda Ind_1 Ind_2 Ind_3 ...];
%     -- for analytic gradient calcs, need to also constrain k to be pos or neg
% 'hold_const' is a list of indices of k that are artificially (and inelegantly held const) for minimization
%
% 'NLtype' specifies spiking non-linearity 0 -> log(1+exp(x)), and 1 -> exp(x).  default is 0

if (nargin < 4) || isempty(silent)
    silent = 0;
end
if (nargin < 5) || isempty(lamrange)
    lamrange = [];
end
if (nargin < 6) || isempty(lamrange2)
    lamrange2 = [];
end
if (nargin < 7) || isempty(Pcon)
    Pcon = [];
end
if (nargin < 8) || isempty(Pmono)
    Pmono = [];
end
if (nargin < 9) || isempty(llist)
    llist = [];
end
if (nargin < 10) || isempty(hold_const)
    hold_const = [];
end
if (nargin < 11) || isempty(NLtype)
    NLtype = 0;
end
if (nargin < 12) || isempty(TolFun)
    TolFun = 1e-7;
end
if nargin < 13
    nlxrange = ones(length(K0)-1,1);
end

[NT NPAR] = size(X);

%if using L1 regularization it must be coupled with constraints on the sign
%of the coefficients
if ~isempty(llist)
    Pcon = [Pcon llist(2:end)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use or construct initial guess %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(K0)
    K0 = zeros(NPAR+1,1);
    K0(1) = 1;
else
    b = K0(end);
    if length(K0) < (NPAR+1)
        K0(end) = 0;
        K0(NPAR+1) = b;
    end
    if length(K0) > (NPAR+1)
        K0 = K0(1:NPAR+1);
        K0(end) = b;
    end
end
K0 = make_row(K0)'; %'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create constraint matrices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Pmono)
    [ND dum] = size(Pmono);
    mults1 = ones(ND,1);
    for n = 1:ND
        if Pmono(n,1) < 0
            mults1(n) = -1;
            Pmono(n,:) = -Pmono(n,:);
        end
    end
    
    Ncon1 = sum(diff(Pmono')); %'
else
    Ncon1 = 0;  ND = 0;
end

Ncon2 = length(Pcon);

if Ncon2 > 0
    mults2 = ones(Ncon2,1);
    mults2(Pcon < 0) = -1;
    Pcon = abs(Pcon);
end

L = length(K0);
A = zeros(Ncon1+Ncon2,L);
B = zeros(Ncon1+Ncon2,1);

ccount = 1;
for n = 1:ND
    for i = Pmono(n,1):(Pmono(n,2)-1)
        A(ccount,[i i+1]) = mults1(n)*[1 -1];
        ccount = ccount + 1;
    end
end

for i = 1:Ncon2
    A(Ncon1+i,Pcon(i)) = -mults2(i);
end

%%%%%%%%%%%%%%%%%%%%%%
% MINIMIZATION SETUP %
%%%%%%%%%%%%%%%%%%%%%%
if silent == 0
    %opts = optimset('GradObj','on','LargeScale','off','Display','iter','MaxIter',2000,'MaxFunEvals',10000,'TolFun',1e-7);
    opts = optimset('GradObj','on','Algorithm','active-set','Display','iter','MaxIter',200,'MaxFunEvals',10000,'TolFun',TolFun);
else
    %opts = optimset('GradObj','on','LargeScale','off','MaxIter',2000,'Display','off','MaxFunEvals',10000,'TolFun',1e-7);
    opts = optimset('GradObj','on','Algorithm','active-set','MaxIter',200,'Display','off','MaxFunEvals',10000,'TolFun',TolFun);
end

%%%% HOLD CONSTANT: things are held constant by making a huge penalty for them to be different
% this is implemented lower-down, so its just a list at this point
if ~isempty(hold_const)
    hold_const = make_row(hold_const);
    hold_const(2,:) = K0(hold_const);
end

% Index of X matrix for spike times: faster to do this here than in LLelog/LLexp ...
SX = X(spksN,:);
Scur_mod_out = cur_mod_out(spksN);

%%%%%%%%%%%%%%%%
% MINIMIZATION %
%%%%%%%%%%%%%%%%
%             [bestk,LL,eflag,outp,grad] = fminunc( @(K) LLelog_jmm(K,X,SX,lamrange,lamrange2,llist,hold_const), K0, opts );
if silent == 0
    options.Display = 'iter';
else
    options.Display = 'off';
end
[bestk,LL,exitflag,output] = minFunc(@(K) LLelog_gradboost(K,X,cur_mod_out,Scur_mod_out,SX,lamrange,lamrange2,llist,hold_const,nlxrange),K0,options);
[LL,grad] = LLelog_gradboost(bestk,X,cur_mod_out,Scur_mod_out,SX,lamrange,lamrange2,llist,hold_const,nlxrange);

%% Calculate final magnitude of gradient as evidence of why fit stopped
grad = sqrt(sum(grad.^2));

%% Make structure with final results
fitp.k = bestk;
fitp.LP = LL;
fitp.LL = LL - LLadjust(bestk,lamrange,length(spksN),lamrange2,llist);
