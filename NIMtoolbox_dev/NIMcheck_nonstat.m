function [fh] = NIMcheck_nonstat(nim,Robs,Xstim,XLin,sm_win)

%%


%% Key parameters
NT = length(Robs);
Nmods = length(nim.mods);
lin_dims = nim.stim_params.lin_dims;
spkhstlen = nim.spk_hist.spkhstlen;
dt = nim.stim_params.dt;

%% make sure Robs is a column vector
if size(Robs,2) > size(Robs,1)
    Robs = Robs';
end

%% CREATE SPIKE HISTORY MATRIX IF NEEDED
if spkhstlen > 0
    Xspkhst = create_spkhist_Xmat(Robs,nim.spk_hist.bin_edges);
else
    Xspkhst = [];
end

%%
theta = nim.spk_NL_params(1); %offset
G = theta + zeros(NT,1); %initialize overall generating function G

Kmat = [nim.mods(:).filtK];
gint = Xstim*Kmat; %subunit generating functions

fgint = nan(NT,Nmods);
for n = 1:Nmods
    
    %process subunit g's with upstream NLs
    if strcmp(nim.mods(n).NLtype,'nonpar')
        fgint(:,n) = piecelin_process(gint(:,n),nim.mods(n).NLy,nim.mods(n).NLx);
    elseif strcmp(nim.mods(n).NLtype,'quad')
        fgint(:,n) = gint(:,n).^2;
    elseif strcmp(nim.mods(n).NLtype,'lin')
        fgint(:,n) = gint(:,n);
    elseif strcmp(nim.mods(n).NLtype,'threshlin')
        fgint(:,n) = gint(:,n);
        fgint(fgint(:,n) < 0,n) = 0;
    else
        error('Invalid internal NL');
    end
    
    %multiply by weight and add to generating function
    G = G + fgint(:,n)*nim.mods(n).sign;
end

%add spike history contribution
if spkhstlen > 0
    G = G + Xspkhst*nim.spk_hist.coefs;
end
%add contribution from linear filter
if lin_dims > 0
    G = G + XLin*nim.kLin;
end
%% Compute predicted firing rate and LL
if strcmp(nim.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    pred_rate = nim.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    pred_rate(too_large) = nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(nim.spk_NL_type,'exp')
    expg = exp(G);
    pred_rate = expg;
else
    error('invalid spk nl');
end
%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(pred_rate) < min_pred_rate
    pred_rate(pred_rate < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%%
tax = (1:length(Robs))*dt;
LLfun = (Robs.* log(pred_rate) - pred_rate); %up to an overall constant
fh = figure;
subplot(2,1,1)
% plot(tax,Robs);
hold on;
plot(tax,smooth(Robs,sm_win),'k','linewidth',1)
plot(tax,smooth(pred_rate,sm_win),'r','linewidth',1)
subplot(2,1,2)
% plot(tax,LLfun);
hold on;
plot(tax,smooth(LLfun,sm_win),'k','linewidth',1)

%%
fprintf('Measured mean rate: %.3f\n',mean(Robs)/dt);
fprintf('Predicted mean rate: %.3f\n',mean(pred_rate)/dt);


