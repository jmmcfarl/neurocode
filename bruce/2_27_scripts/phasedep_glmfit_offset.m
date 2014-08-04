function [LL, LLgrad] = phasedep_glmfit_offset(params, fix_params, Robs, Xstim,XPBx,XLin,L2_mat,reg_lambda,stim_params,PBx)
%

%% USEFUL PARAMETERS
filtLen = prod(stim_params.stim_dims);
nBars = squeeze(stim_params.stim_dims(2));
nPBs = length(PBx);
lin_dims = stim_params.lin_dims;

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
% theta = params(end); %offset
kmat = reshape(fix_params,filtLen,nPBs);
KPLin = params(1:nPBs);
KLin = params((nPBs+1):(nPBs+lin_dims));
filt_outs = Xstim*kmat;

% G = XPBx*Klin + theta;%initialize overall generating function G
G = XPBx*KPLin + XLin*KLin;%initialize overall generating function G

%%
G = G + sum(filt_outs.*XPBx,2);
%% Compute predicted firing rate
max_gbeta = 50; %to prevent numerical overflow
expg = exp(G);
too_large = (G > max_gbeta);
r = log(1+expg); %alpha*log(1+exp(gbeta))
r(too_large) = G(too_large); %log(1+exp(x)) ~ x in limit of large x

%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(r) < min_pred_rate
    r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%% COMPUTE LL and LL gradient
LL = sum(Robs.* log(r) - r); %up to an overall constant

%'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
residual = (Robs./r - 1) .* expg ./ (1+expg);
residual(too_large) = (Robs(too_large)./r(too_large) - 1);

%initialize LL gradient
LLgrad = zeros(length(params),1)';

% % Calculate derivatives with respect to constant term (theta)
% LLgrad(end) = sum(residual);

LLgrad((1):(nPBs)) = residual'*XPBx;
LLgrad((nPBs + 1):(nPBs+lin_dims)) = residual'*XLin;

%%
%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad'/Nspks;

