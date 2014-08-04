function fo = fitWeights_stcb_nonlpsc(fo,g_mat,spkbs,silent,TolFun)
%[fo, LL] = fitWeights(fo,input,spkbs)
% fitWeights: fit modul-weights and intercept of the full GNLM

if nargin < 4
    silent = 1;
end
if nargin < 5
    TolFun = 1e-7;
end

nmods = length(fo.mods); 

rect_gmat = g_mat;
rect_gmat(rect_gmat < 0) = 0;

fX = [];
for i = 1:nmods
%     fX = [fX nlin_proc_stim(g_mat(:,i),fo.mods(i).nly,fo.mods(i).nlx)];
    fX = [fX rect_gmat(:,i)];
end

%compile all module weights into a single vector
W0 = [];
for i=1:nmods
    W0=[W0; fo.mods(i).w];
end
W0 = [W0; fo.const]; %tack on the constant term
 
%initialize parameters
lamrange = [];
lamrange2 = [];
Pcon = [];
Pmono = [];
hold_const = [];
NLtype = 0;

% %for L1 regularization of weights
% llist = [];
% if fo.lambdaW > 0
%     llist = [fo.lambdaW];
%     llist = [llist 1:nmods]; %all constrained to be positive
% end
% [fitp,grad] = GLMsolve_jmm( fX, spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype,TolFun );
% fo.const = fitp.k(end);
% for imod = 1:nmods
%     fo.mods(imod).w = fitp.k(imod);
% end
% % fo.LL = fitp.LL; 
% % fo.LLnadj = fitp.LLnadj;
% end

% %compute binned spike vector
% Robs = zeros(1,size(g_mat,1));
% ftable = tabulate(spkbs);
% Robs(ftable(:,1)) = ftable(:,2);
lambda = ones(1,nmods)*fo.lambdaW;
lambda = [lambda 0]/length(spkbs);
SX = fX(spkbs,:);
% if max(lambda) > 0

if strcmp(fo.spk_nl,'logexp')
fitp = L1General2_PSSas(@(K) LLelog_jmm(K,fX,SX,[],[],[],[]),W0,lambda',[],1);
elseif strcmp(fo.spk_nl,'exp')
  fitp = L1General2_PSSas(@(K) LLeexp_jmm(K,fX,SX,[],[],[],[]),W0,lambda',[],1);  
end
fo.const = fitp(end);
for imod = 1:nmods
    fo.mods(imod).w = fitp(imod);
end



