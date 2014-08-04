function [fo] = fitHI_jmm(fo,input,spkbs,silent)
% fitHI: fit history term for internal receptive field

flen  = length(fo.mods(1).h); %length of PSC term 
nmods = length(fo.mods); %number of modules


%% apply nonlinearity to the filtered stimulus -> f(kx)
fgs = zeros(size(input)); %initialize f(kx) matrix 
for imod= 1:nmods; 
    fgs(:,imod) = nlin_proc_stim(input(:,imod), fo.mods(imod).nly, fo.mods(imod).nlx);
end;

NNT = size(fgs,1); 
%% create time embedding for f(kx)
%is there any way to do this where the loop is over 1:flen rather than over
%all time samples??
X  = zeros(NNT-flen+1,nmods*flen);
for i = flen:NNT
    for j =1:nmods 
    X(i-flen+1,((j-1)*flen+1):(j*flen)) = fgs((i-flen+1):i,j);
    end
end

% S  = X(spkbs,:); %X evaluated at spike times

%compile all history terms into a single vector across modules
IH0 = []; 
for i=1:nmods
%     IH0=[IH0; fo.mods(i).w*fo.mods(i).h];
      IH0=[IH0; fo.mods(i).h]; %do not incorporate module weighting into K vector for regularization purposes
end 
IH0 = [IH0; fo.const]; %add in constant term

%% compile lambda range matrix
lamrange = [];
for imod=1:nmods
    if fo.mods(imod).lh ~= 0
        lamrange = [lamrange; fo.mods(imod).lh (imod-1)*flen+1 imod*flen];
    end
end

%% create lamrange matrix
lamrange2 =[];
for i=1:nmods
    if fo.mods(i).lh2 ~= 0
        lamrange2 = [lamrange2; fo.mods(i).lh2, (i-1)*flen+1, i*flen];
    end
end
%% create llist matrix (l1 penalty)
llist = [];
% for i = 1:nmods
%     if fo.mods(i).lh1 ~= 0
%         llist = [llist ((i-1)*flen+1):(i*flen)];
%         cur_lh1 = fo.mods(i).lh1;
%     end
% end
% if ~isempty(llist)
%     llist = [cur_lh1 llist];
% end

%%
% silent = 0;
Pcon = [];
for i = 1:nmods
    if fo.mods(i).hcon ~= 0
        Pcon = [Pcon fo.mods(i).hcon*(((i-1)*flen+1):i*flen)];
    end
end
Pmono = [];
for i = 1:nmods
    if fo.mods(i).hmon ~= 0
        Pmono = [Pmono; fo.mods(i).hmon*[((i-1)*flen+1) i*flen]];
    end
end
hold_const = [];
NLtype = 0;
[fitp,grad] = GLMsolve_jmm( X, spkbs, IH0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype );
fo.const = fitp.k(end);
himat = reshape(fitp.k(1:end-1),flen,nmods); %reshape best K vector into a matrix of PSC coefs

fo.LL = fitp.LL; 
fo.LP = fitp.LP;

%% normalize history term
for i=1:nmods
    thish        = himat(:,i);
    hsum         = sum(abs(thish));
    fo.mods(i).h = thish/hsum; %normalize h coefs to have sum 1
    fo.mods(i).w = hsum;
%     fo.mods(i).w = fo.mods(i).w*hsum; %pull the normalization factor INTO the weight term
end

end

