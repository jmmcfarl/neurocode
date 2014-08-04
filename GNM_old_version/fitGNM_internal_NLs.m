function [fo] = fitGNM_internal_NLs(fo,Stim,spkbs,silent,adj)

if nargin < 5
    adj = 1; %default is that wew will adjust the scale of the internal NL
end

stimlen   = size(Stim,1);  %number of samples of filtered stim
g = Stim*get_k_mat(fo);
nmods = length(fo.mods); %number of modules
NBFS  = length(fo.mods(1).nlx); %number of basis functions for non-linearity

X = zeros(stimlen,nmods*NBFS); %initialize X matrix which is for the NL BFs of each module
for imod=1:nmods; 	%for each module
    
%     %update the tent basis X-axis
%     n_dxs = length(fo.mods(imod).nlx);
% %         fo.mods(imod).nlx = prctile(g(:,imod),linspace(0.05,99.95,n_dxs)); %equipopulated
% %     fo.mods(imod).nlx = linspace(prctile(g(:,imod),0.05),prctile(g(:,imod),99.95),n_dxs); %equispacing
%     
%     %set nearest tent basis to 0
%     [~,nearest] = min(abs(fo.mods(imod).nlx));
%     fo.mods(imod).nlx(nearest) = 0;
%     fo.mods(imod).nly = fo.mods(imod).nlx;
%     fo.mods(imod).nly(1:nearest) = 0;
    
    % the output of the current model's internal filter projected onto the tent basis representation
    tbfi = fo.mods(imod).w*tbrep(g(:,imod),fo.mods(imod).nlx);
    
    X(: , ((imod-1)*NBFS + 1) : (imod*NBFS)) = tbfi; %assemble filtered NLBF outputs into X matrix
end;

%%
% create lamrange matrix
lamrange =[];
for i=1:nmods
    lamrange = [lamrange; fo.mods(i).lnl, (i-1)*NBFS+1, i*NBFS];
end

% create lamrange matrix
lamrange2 =[];
for i=1:nmods
    lamrange2 = [lamrange2; fo.mods(i).lnl2, (i-1)*NBFS+1, i*NBFS];
end

% compile all NLBF coefs into a single vector
nly = [];
mod_inds = [];
for i=1:nmods
    nly = [nly fo.mods(i).nly]; %not incorporating the multiplier here because doing so messes with regularization
    mod_inds = [mod_inds i*ones(size(fo.mods(i).nly))];
end
% nly = [nly fo.const]'; %tack on the constant term
nly = [nly fo.spk_theta]'; %tack on the constant term

Pcon = [];
for i = 1:nmods
    if fo.mods(i).nlcon ~= 0
        Pcon = [Pcon fo.mods(i).nlcon*(((i-1)*NBFS+1):i*NBFS)];
    end
end

%contstrain NL coefs to be monotonically increasing for each module
Pmono = [];
for i = 1:nmods
    if fo.mods(i).nlmon ~= 0
        Pmono = [Pmono; fo.mods(i).nlmon*[((i-1)*NBFS+1) i*NBFS]];
    end
end
llist = [];

%constrain the NL to be 0 at 0 input
Pzero = [];
n_pts = 0;
for i = 1:nmods
    zp = find(fo.mods(i).nlx==0);
    Pzero = [Pzero zp+n_pts];
    n_pts = n_pts + length(fo.mods(i).nlx);
end

hold_const = [];
NLtype = 0; %log(1+exp(g))
% [fitp,grad] = GLMsolve_jmm( X, spkbs, nly, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype ,1e-5,[],Pzero);
[fitp,grad] = GLMsolve_gnm(fo, X, spkbs, nly, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const ,1e-5,[],Pzero);

%%
nlmat    = reshape(fitp.k(1:end-1),NBFS,nmods); %take output K vector and restructure into a matrix of NLBF coefs, one for each module
nlmat_fin = nlmat;
nlx = fo.mods(1).nlx';
for i=1:nmods;
    cur_pset = ((i-1)*NBFS+1) : (i*NBFS);
    thisnl = nlmat(:,i); %NL coefs for current module
    
    %if you want to rescale the internal NL coefs to match the scale before
    %fitting
    if adj == 1
      prev_range = range(fo.mods(i).nly);
      %         thisnl = thisnl-min(thisnl);
        nlrange = range(thisnl);
        thisnl = thisnl/nlrange*prev_range;
    elseif adj == 2
        prev_std = std(X(:,cur_pset)*fo.mods(i).nly');
        cur_std = std(X(:,cur_pset)*thisnl);
        thisnl = thisnl*prev_std/cur_std;
    end
    fo.mods(i).nly = thisnl';
    nlmat_fin(:,i) = fo.mods(i).nly';
end

fin_nlys = nlmat_fin(:);
fin_nlys(isnan(fin_nlys) | isinf(fin_nlys)) = 0;
g = X*fin_nlys;
%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spkbs);
Robs(ftable(:,1)) = ftable(:,2);

gs = g(spkbs);
if silent == 0
    opts = optimset('GradObj','on','Display','iter','MaxIter',500,'MaxFunEvals',10000,'TolFun',1e-10);
else
    opts = optimset('GradObj','on','MaxIter',500,'Display','off','MaxFunEvals',10000,'TolFun',1e-10);
end
[best_const,LL,eflag] = fminunc( @(K) internal_theta_opt(K,g,Robs,fo), mean(g), opts );
fo.spk_theta = best_const;

fo.LL = fitp.LL;
fo.LP = fitp.LP;
end

function [ll,grad] = internal_theta_opt(b,g,Robs,fo)

eg = exp(fo.spk_beta*(g+b));
r = fo.spk_alpha*log(1+eg);
r(r < 1e-20) = 1e-20;
% Standard point-process likelihood

ll = sum(Robs'.*log(r) - r);

residual = fo.spk_alpha * fo.spk_beta * (Robs'./r - 1) .* (eg./(eg + 1));
grad = sum(residual);

% wg = eg./(1+eg);
% ers = (egs./(1+egs))./rs;
% % grad = ers'*gs-wg'*g;             
% grad = sum(ers)-sum(wg);  % grad w.r.t. b

nspks = sum(Robs); %normalize by number of spikes
ll=-ll/nspks;
grad=-grad'/nspks; %'

end

