function [fitp,grad,se] = GLMsolve_jmm( X, spksN, K0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, TolFun, nlxrange, Pzero)
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
if nargin < 13 || isempty(nlxrange)
    nlxrange = ones(length(K0),1);
end
if nargin < 14
    Pzero = [];
end

[NT NPAR] = size(X);

% %if using L1 regularization it must be coupled with constraints on the sign
% %of the coefficients
% if ~isempty(llist)
%     Pcon = [Pcon llist(2:end)];
% end

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
Ncon3 = length(Pzero);

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

Aeq = zeros(Ncon3,L);
Beq = zeros(Ncon3,1);
for i = 1:Ncon3
    Aeq(i,Pzero(i)) = 1;
end


%%%%%%%%%%%%%%%%%%%%%%
% MINIMIZATION SETUP %
%%%%%%%%%%%%%%%%%%%%%%
if silent == 0
    opts = optimset('GradObj','on','LargeScale','off','Display','iter','MaxIter',2000,'MaxFunEvals',10000,'TolFun',1e-7);
    %     opts = optimset('GradObj','on','Algorithm','active-set','Display','iter','MaxIter',200,'MaxFunEvals',10000,'TolFun',TolFun);
else
    opts = optimset('GradObj','on','LargeScale','off','MaxIter',2000,'Display','off','MaxFunEvals',10000,'TolFun',1e-7);
    %     opts = optimset('GradObj','on','Algorithm','active-set','MaxIter',200,'Display','off','MaxFunEvals',10000,'TolFun',TolFun);
end

%%%% HOLD CONSTANT: things are held constant by making a huge penalty for them to be different
% this is implemented lower-down, so its just a list at this point
if ~isempty(hold_const)
    hold_const = make_row(hold_const);
    hold_const(2,:) = K0(hold_const);
end

% Index of X matrix for spike times: faster to do this here than in LLelog/LLexp ...
SX = X(spksN,:);

%%%%%%%%%%%%%%%%
% MINIMIZATION %
%%%%%%%%%%%%%%%%
if NLtype == 0
    if (Ncon1 + Ncon2 + Ncon3) > 0
        [bestk,LL,eflag,outp,lam,grad] = fmincon( @(K) LLelog_jmm(K,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange), K0, A, B, Aeq, Beq,[],[],[], opts );
    else
        %             [bestk,LL,eflag,outp,grad] = fminunc( @(K) LLelog_jmm(K,X,SX,lamrange,lamrange2,llist,hold_const), K0, opts );
        if silent == 0
            options.Display = 'iter';
        else
            options.Display = 'off';
        end
        if ~isempty(llist)
            
            lambda = zeros(size(K0));
            lambda(llist(2:end)) = llist(1);
            lambda = lambda/length(spksN);
            [bestk LL] = L1General2_PSSas(@(K) LLelog_jmm(K,X,SX,lamrange,lamrange2,[],hold_const,nlxrange),K0,lambda,options);
        else
            options.progTol = 1e-10;
            [bestk,LL,exitflag,output] = minFunc(@(K) LLelog_jmm(K,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange),K0,options);
        end
        [LL,grad] = LLelog_jmm(bestk,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange);
        
        if nargout > 2
            disp('Computing standard errors')
            H = nan(length(bestk));
            h = 1e8*eps;
            for i = 1:length(bestk)
                % derivative at first point (left)
                x1 = bestk;
                x1(i) = x1(i) - h;
                [LL,grad1] = LLelog_jmm(x1,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange);
                
                % derivative at second point (right)
                x2 = bestk;
                x2(i) = x2(i) + h;
                [LL,grad2] = LLelog_jmm(x2,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange);
                
                % differentiate between the two derivatives
                d2f = (grad2-grad1)*length(spksN) / (2*h);
                
                % assign as row i of Hessian
                H(i,:) = d2f';
                
            end
            invH = pinv(H);
            se = sqrt(diag(invH));
        end
    end
else
    %%% Exponential non-linearity
    if (Ncon1 + Ncon2 + Ncon3) > 0
        [bestk,LL,eflag,outp,lam,grad] = fmincon( @(K) LLexp_jmm(K,X,SX,lamrange,lamrange2,llist,hold_const), K0, A, B, [],[],[],[],[], opts );
    else
        if silent == 0
            options.Display = 'iter';
        else
            options.Display = 'off';
        end
        if ~isempty(llist)
            lambda = zeros(size(K0));
            lambda(llist(2:end)) = llist(1);
            lambda = lambda/length(spksN);
            [bestk LL] = L1General2_PSSas(@(K) LLexp_jmm(K,X,SX,lamrange,lamrange2,[],hold_const,nlxrange),K0,lambda,options);
        else
            %       [bestk,LL,eflag,outp,grad] = fminunc( @(K) LLexp_jmm(K,X,SX,lamrange,lamrange2,llist,hold_const), K0, opts );
            [bestk,LL,exitflag,output] = minFunc(@(K) LLexp_jmm(K,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange),K0,options);
        end
        [LL,grad] = LLexp_jmm(bestk,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange);
        
        if nargout > 2
            disp('Computing standard errors')
            H = nan(length(bestk));
            h = 1e3*eps;
            for i = 1:length(bestk)
                % derivative at first point (left)
                x1 = bestk;
                x1(i) = x1(i) - h;
                [LL,grad1] = LLexp_jmm(x1,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange);
                
                % derivative at second point (right)
                x2 = bestk;
                x2(i) = x2(i) + h;
                [LL,grad2] = LLexp_jmm(x2,X,SX,lamrange,lamrange2,llist,hold_const,nlxrange);
                
                % differentiate between the two derivatives
                d2f = (grad2-grad1)*length(spksN) / (2*h);
                
                % assign as row i of Hessian
                H(i,:) = d2f';
                
            end
            invH = inv(H);
            se = sqrt(diag(invH));
        end
    end
end
%% Calculate final magnitude of gradient as evidence of why fit stopped
grad = sqrt(sum(grad.^2));

%% Make structure with final results
fitp.k = bestk;
fitp.LP = LL;
fitp.LL = LL - LLadjust(bestk,lamrange,length(spksN),lamrange2,llist);
