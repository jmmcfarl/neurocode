function [fo] = fitNL_jmm(fo,input,spkbs,silent)
% USAGe: [fo] = fitNL(fo,input,spkbs)
%fitNL: fit nonllinearity by adjusting bf-coefficients
% fo:     the fittable object,
% input:  the binned input,
% spkbs:  timebins where spikes occurred

%input is the stimulus filtered by each modules RF

hilen = length(fo.mods(1).h); %number of lags in the PSC term
NNT   = size(input,1);  %number of samples of filtered stim
nmods = length(fo.mods); %number of modules
NBFS  = length(fo.mods(1).nlx); %number of basis functions for non-linearity

X = zeros(NNT-hilen+1,nmods*NBFS); %initialize X matrix which is for the NL BFs of each module

for imod=1:nmods; 	%for each module
    
    % the output of the current model's internal filter projected onto the tent basis representation
    tbfi = tbrep(input(:,imod),fo.mods(imod).nlx);
    
    % apply temporal filter
    ftb = zeros(NNT-hilen+1,NBFS); %initialize output of temporal filter
    %for each NL basis function, convolve its output with the PSC kernel
    for ibf = 1:NBFS
        % 	  fgh = filter(fo.mods(imod).h(end:-1:1),1,tbfi(:,ibf)); %filter rather than conv?
        % 	  ftb(:,ibf) = fgh(hilen:end);
        
        %SHOULD WE INCLUDE WEIGHTINGS HERE OR NOT??
%         fgh = g_convolve_jmm(tbfi(:,ibf),fo.mods(imod).h(end:-1:1));
               fgh = fo.mods(imod).w*g_convolve_jmm(tbfi(:,ibf),fo.mods(imod).h(end:-1:1));
        ftb(:,ibf) = fgh(hilen:end);
    end
    
    X(: , ((imod-1)*NBFS + 1) : (imod*NBFS)) = ftb; %assemble filtered NLBF outputs into X matrix
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
for(i=1:nmods)
    %     nly = [nly,fo.mods(i).w*fo.mods(i).nly];
    nly = [nly fo.mods(i).nly]; %not incorporating the multiplier here because doing so messes with regularization
end
nly = [nly fo.const]'; %tack on the constant term
% silent = 0;
Pcon = [];
for i = 1:nmods
    if fo.mods(i).nlcon ~= 0
        Pcon = [Pcon; fo.mods(i).nlcon*(((i-1)*NBFS+1):i*NBFS)];
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
hold_const = [];
NLtype = 0;
[fitp,grad] = GLMsolve_jmm( X, spkbs, nly, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype );

%%
fo.const = fitp.k(end);
nlmat    = reshape(fitp.k(1:end-1),NBFS,nmods); %take output K vector and restructure into a matrix of NLBF coefs, one for each module

for i=1:nmods;
    thisnl         = nlmat(:,i); %NL coefs for current module
    nlorigin       = find(fo.mods(i).nlx==0); %find origin of NL x-ax
    
    %shift coef values so that f(0)=0
%     offset         = thisnl(nlorigin);
%     shifted        = thisnl - offset;
%     
%     nlrange        = range(shifted); %compute the range of coef values
%     fo.mods(i).nly = (shifted/nlrange)'; %scale so the range is 1 WHY??
fo.mods(i).nly = thisnl;

    %     fo.mods(i).w = fo.mods(i).w*nlrange;
    %     fo.mods(i).nly = shifted';
    % weights are fit by fitHI...
end;
fo.LL = fitp.LL;
fo.LP = fitp.LP;

