clear all;close all


%% CREATE A 1d binary noise stimulus like in Rust et al., 2005
dt = 0.01; %in s
NT = 50000; %number time samples
xvNT = 50000; %number time samples for simulated cross-val
SDIM = 24; %spatial dimensions
flen = 14; %time lags
stim = round(2*rand(NT,SDIM))-1;
xvstim = round(2*rand(xvNT,SDIM))-1;

%% GENERATE FILTERS
nfilts = 4;
x = repmat(1:SDIM,flen,1); %x-axis 

LAMBDA = 5; %spatial freq
SIGMA = 1.75; %spatial std dev
sigma1 = repmat(SIGMA,flen,SDIM); %constant spatial width
lambda = repmat(LAMBDA,flen,SDIM); %constant SF
amp_vec = [0 0 gampdf((1:flen-2)-1,3,1.3)]; %temporal weighting function
amp_vec = fliplr(amp_vec); 
weight_vec = repmat(amp_vec',1,SDIM);
offset_vec = zeros(flen,SDIM); %no offset

spacing = LAMBDA; %spacing between center positions
xovals = 0:spacing:(spacing*(nfilts-1));
xovals = xovals - mean(xovals);
xovals = xovals + SDIM/2+0.5;

filt_mat = zeros(nfilts,flen,SDIM);
filt_vec = zeros(nfilts,flen*SDIM);
for i = 1:nfilts
    psi1 = repmat(((1:flen)'/6).^3+5,1,SDIM); %make direction selective by having the spatial phase evolve in time
    
    x0 = repmat(xovals(i),flen,SDIM);
    temp = weight_vec.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+offset_vec;
    filt_mat(i,:,:) = temp/norm(temp(:));
    filt_vec(i,:) = temp(:)/norm(temp(:));
end
filt_vec = filt_vec';


%% CREATE TIME-EMBEDDED STIMULUS
stim_emb = makeStimRows(stim,flen);
xvstim_emb = makeStimRows(xvstim,flen);
% FILTER STIMULUS
filt_stim = zeros(nfilts,NT);
xvfilt_stim = zeros(nfilts,xvNT);
for i = 1:nfilts
    filt_stim(i,:) = zscore(stim_emb*filt_vec(:,i));
    xvfilt_stim(i,:) = zscore(xvstim_emb*filt_vec(:,i));
end


%% create cond intens function
%pass through log(1+exp(x)) internal NLs
beta = 4; theta = 0.5; %internal NL parameters
g = zeros(1,NT); %gen function
xvg = zeros(1,xvNT);

%PASS OUTPUT OF EACH FILTER THROUGH INTERNAL NLs
for i = 1:nfilts
    lfilt_stim(i,:) = 1/beta*log(1+exp(beta*(filt_stim(i,:)-theta)));
    xvlfilt_stim(i,:) = 1/beta*log(1+exp(beta*(xvfilt_stim(i,:)-theta)));
    g = g + lfilt_stim(i,:);
    xvg = xvg + xvlfilt_stim(i,:);
end
%normalize outputs
sg = std(g);
g = g/sg;
xvg = xvg/sg;

%APPLY SPIKING NL (log(1+exp(x)))
cur_theta = 1; 
cur_beta = 4;
p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
xvp_spike = 1/cur_beta*log(1+exp(cur_beta*(xvg-cur_theta)));
target_rate = 20; %in Hz
target_pspike = target_rate*dt;
scale_f = target_rate/(mean(p_spike)/dt);
p_spike = p_spike*scale_f;
xvp_spike = xvp_spike*scale_f;

%% GENERATE SPIKES
spikes = poissrnd(p_spike);
spikebins = convert_to_spikebins(spikes);
fprintf('Nspks: %d\n',length(spikebins));

xvspikes = poissrnd(xvp_spike);
xvspikebins = convert_to_spikebins(xvspikes);

avg_rate = length(spikebins)/NT;
prate = ones(1,xvNT)*avg_rate;
null_xvLL = -sum(xvspikes.*log(prate)-prate)/length(xvspikebins); %cross-validated likelihood of the null model

%% QUAD MOD
%structure containing stimulus parameters
disp('Fitting quadratic model')
stim_params.spatial_dims = 1; %dimensionality of stimulus
stim_params.sdim = SDIM; %number of spatial pixels
stim_params.flen = flen;

defmod.lambda_L1x = 5; %L1 penalty on filter coeffs
defmod.lambda_d2XT = 5; %Spatiotemporal laplacian penalty (equal in spatial and temporal dims)

nmods = 4; %number of subunits
%one linear and the rest squared nonlinearities
kern_types{1} = 'lin';
for j = 2:nmods
    kern_types{j} = 'quad';
end
%random intialization on filter coefs
init_kerns = randn(flen*SDIM,nmods);
init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2))); %rescale

%set which filters are 'excitatory' (+1) and 'suppressive' (-1)
init_signs = ones(1,nmods); %all exc.
quad_mod = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params); %initialize model

%estimate filter coefs
quad_mod = fitGNM_filters(quad_mod,stim_emb,spikebins,'none',[],1e-4,1e-6);

%estimate spk NL params
[~, ~, ~, pred_rate, g] = getLL_GNM(quad_mod,stim_emb,spikebins,'none');
quad_mod = fitGNM_spkNL(quad_mod,g,spikebins,0);

%cross-val log-like
quad_xvLL = getLL_GNM(quad_mod,xvstim_emb,xvspikebins,'none')

%% GNM
disp('Fitting GNM')
nmods = 4; %number of subunits
%one linear and the rest squared nonlinearities
for j = 1:nmods
    kern_types{j} = 'threshlin'; %initialize internal NLs as threshold linear functions
end
defmod.lambda_L1x = 100; %L1 penalty on filter coeffs
defmod.lambda_d2XT = 10; %Spatiotemporal laplacian penalty (equal in spatial and temporal dims)
defmod.lnl2 = 20; %L2 penalty on smoothness of upstream NL coefs
defmod.nlmon = 1; %constrain upstream NLs to be monotonic


%random intialization on filter coefs
init_kerns = randn(flen*SDIM,nmods); init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
init_signs = ones(1,nmods);
gnm = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);

%fit filters
gnm = fitGNM_filters(gnm,stim_emb,spikebins,'none',[],1e-4,1e-6);

%switch to using tent-basis representation of upstream NLs
gnmr = setGNM_NLBFs(gnm,stim_emb); %this sets the tent-basis functions based on the generating distributions
gnmr = adjust_all_reg(gnmr,'nltype','uncon');

%fit upstream NLs
gnmr = fitGNM_internal_NLs(gnmr,stim_emb,spikebins,1,2);

%fit spiking NL
[~, ~, ~, ~, g] = getLL_GNM(gnmr,stim_emb,spikebins,'none');
gnmr = fitGNM_spkNL(gnmr,g,spikebins,0);

%re-estimate filters
gnmr = fitGNM_filters(gnmr,stim_emb,spikebins,'none',[],1e-4,1e-6);

%re-estimate upstream NLs
gnmr = fitGNM_internal_NLs(gnmr,stim_emb,spikebins,1,2);

gnmr_xvLL = getLL_GNM(gnmr,xvstim_emb,xvspikebins,'none')

%% VISUALIZE MODEL FITS
n_rows = 4;
sc_type = 'centered';
plotfo1d_nopsc(quad_mod,n_rows,sc_type)
plotfo1d_nopsc(gnmr,n_rows,sc_type)
