function [fo] = fitNL_lexp_adj(fo,input,spkbs,silent,adj)
% USAGe: [fo] = fitNL(fo,input,spkbs)
%fitNL: fit nonllinearity by adjusting bf-coefficients
% fo:     the fittable object,
% input:  the binned input,
% spkbs:  timebins where spikes occurred

%input is the stimulus filtered by each modules RF

if nargin < 5
    adj = 1;
end

hilen = 1; %number of lags in the PSC term
NNT   = size(input,1);  %number of samples of filtered stim
nmods = length(fo.mods); %number of modules
NBFS  = length(fo.mods(1).nlx); %number of basis functions for non-linearity

X = zeros(NNT,nmods*NBFS); %initialize X matrix which is for the NL BFs of each module

for imod=1:nmods; 	%for each module
    
    %update the tent basis X-axis
    n_dxs = length(fo.mods(imod).nlx);
    %     fo.mods(imod).nlx = prctile(input(:,imod),linspace(1,99,n_dxs)); %equipopulated
    fo.mods(imod).nlx = linspace(prctile(input(:,imod),0.25),prctile(input(:,imod),99.75),n_dxs); %equispacing
    
    %set nearest tent basis to 0
    [~,nearest] = min(abs(fo.mods(imod).nlx));
    fo.mods(imod).nlx(nearest) = 0;
    
    % the output of the current model's internal filter projected onto the tent basis representation
    tbfi = fo.mods(imod).w*tbrep(input(:,imod),fo.mods(imod).nlx);
    %     tbfi = tbrep(input(:,imod),fo.mods(imod).nlx);
    
    X(: , ((imod-1)*NBFS + 1) : (imod*NBFS)) = tbfi; %assemble filtered NLBF outputs into X matrix
end;

% S   = X(spkbs,:); %precompute X-matrix at spikes

%% create lamrange matrix
lamrange =[];
for i=1:nmods
    lamrange = [lamrange; fo.mods(i).lnl, (i-1)*NBFS+1, i*NBFS];
end

%% create lamrange matrix
lamrange2 =[];
for i=1:nmods
    lamrange2 = [lamrange2; fo.mods(i).lnl2, (i-1)*NBFS+1, i*NBFS];
end

%% solve for updated params using GLMsolve_jmm
%compile all NLBF coefs into a single vector
nly = [];
mod_inds = [];
for(i=1:nmods)
    %     nly = [nly,fo.mods(i).w*fo.mods(i).nly];
    nly = [nly fo.mods(i).nly]; %not incorporating the multiplier here because doing so messes with regularization
    mod_inds = [mod_inds i*ones(size(fo.mods(i).nly))];
end
nly = [nly fo.const]'; %tack on the constant term
% silent = 0;
Pcon = [];
for i = 1:nmods
    %     fo.mods(i).nlcon = 1;
    if fo.mods(i).nlcon ~= 0
        Pcon = [Pcon fo.mods(i).nlcon*(((i-1)*NBFS+1):i*NBFS)];
    end
end

%contstrain NL coefs to be monotonically increasing for each module
Pmono = [];
for i = 1:nmods
    %     fo.mods(i).nlmon = 1;
    if fo.mods(i).nlmon ~= 0
        Pmono = [Pmono; fo.mods(i).nlmon*[((i-1)*NBFS+1) i*NBFS]];
    end
end
llist = [];

Pzero = [];
n_pts = 0;
for i = 1:nmods
    zp = find(fo.mods(i).nlx==0);
    Pzero = [Pzero zp+n_pts];
    n_pts = n_pts + length(fo.mods(i).nlx);
end

%create list of parameters which are held constant
hold_const = [];
% for i = 1:nmods
%    if ~strcmp(fo.mods(i).nltype,'uncon')
%        hold_const = [hold_const find(mod_inds==i)];
%    end
% end

NLtype = 0;
[fitp,grad] = GLMsolve_jmm( X, spkbs, nly, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype ,1e-5,[],Pzero);

%%
% fo.const = fitp.k(end);
nlmat    = reshape(fitp.k(1:end-1),NBFS,nmods); %take output K vector and restructure into a matrix of NLBF coefs, one for each module
nlmat_fin = nlmat;
nlx = fo.mods(1).nlx';
for i=1:nmods;
    prev_range = range(fo.mods(i).nly);
    cur_pset = ((i-1)*NBFS+1) : (i*NBFS);
    prev_scale = std(X(:,cur_pset)*nly(cur_pset));
    
    thisnl         = nlmat(:,i); %NL coefs for current module
    
    if adj == 1
%         thisnl = thisnl-min(thisnl);
            nlrange = range(thisnl);
            thisnl = thisnl/nlrange*prev_range;
        
%         new_scale = std(X(:,cur_pset)*thisnl);
%         thisnl = thisnl/new_scale*prev_scale;
        
     end
       fo.mods(i).nly = thisnl';
        %     fo.mods(i).w = fo.mods(i).w*nlrange; %pull NL range into w
        %     fit_beta = lsqnonlin(@(b) logexp(nlx,b)-thisnl,1,0.01,50);
        %     fo.mods(i).beta = fit_beta;
        %     fo.mods(i).nly = logexp(nlx,fit_beta)';
    nlmat_fin(:,i) = fo.mods(i).nly';
    
end;


fin_nlys = nlmat_fin(:);
g = X*fin_nlys;
gs = g(spkbs);
if silent == 0
    %   opts = optimset('GradObj','on','Algorithm','active-set','Display','iter','MaxIter',500,'MaxFunEvals',10000,'TolFun',1e-10);
    opts = optimset('GradObj','on','Display','iter','MaxIter',500,'MaxFunEvals',10000,'TolFun',1e-10);
else
    %   opts = optimset('GradObj','on','Algorithm','active-set','MaxIter',500,'Display','off','MaxFunEvals',10000,'TolFun',1e-10);
    opts = optimset('GradObj','on','MaxIter',500,'Display','off','MaxFunEvals',10000,'TolFun',1e-10);
end
[best_const,LL,eflag] = fminunc( @(K) const_ll_opt(K,g,gs), mean(g), opts );
fo.const = best_const;

fo.LL = fitp.LL;
fo.LP = fitp.LP;

end
