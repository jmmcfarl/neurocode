cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
% close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

NT = 100000; xvNT = 50000; 
SDIM = 24; flen = 14;
% stim = round(2*rand(NT,SDIM))-1;
% xvstim = round(2*rand(xvNT,SDIM))-1;
stim = trnd(10,NT,SDIM);
xvstim = trnd(10,xvNT,SDIM);

% GENERATE A FILTER
nfilts = 6;
LAMBDA = 5; %5 4
SIGMA = 1.5; %1.5
x = repmat(1:SDIM,flen,1);
sigma1 = repmat(SIGMA,flen,SDIM);
lambda = repmat(LAMBDA,flen,SDIM);
% amp_vec = ncx2pdf((1:flen-2)-1,4,1);
% amp_vec = ncx2pdf((1:flen-2)-1,4,0);
amp_vec = gampdf((1:flen-2)-1,3,1.3);
amp_vec = [0 0 amp_vec];
amp_vec = fliplr(amp_vec);
b = repmat(amp_vec',1,SDIM);
a = repmat(0,flen,SDIM);

desired_spacing = LAMBDA/5.75; %5.5/ 5 / 4.5n
beg = SDIM/2-desired_spacing*floor(nfilts/2)+0.75;
xovals = 0:desired_spacing:(desired_spacing*(nfilts-1));
xovals = xovals - mean(xovals);
xovals = xovals + SDIM/2+0.5;
for i = 1:nfilts
    psi1 = repmat(((1:flen)'/6).^3+5,1,SDIM); 

    x0 = repmat(xovals(i),flen,SDIM);
%     cur_psi = mod(psi1 + rand*pi,2*pi);
    temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
    filt(i,:,:) = temp/norm(temp(:));
    filt_mat(i,:) = temp(:)/norm(temp(:));
end
filt_mat = filt_mat';
c_cent = SDIM/2+0.5;
c_std = 1.6; %1.6
max_cval = 0.5; %1
c_offset = 1;
cvals = max_cval*exp(-(xovals-c_cent).^2/(2*c_std))+c_offset;

% CREATE TIME-EMBEDDED STIMULUS
stim_emb = makeStimRows(stim,flen);
xvstim_emb = makeStimRows(xvstim,flen);

% FILTER STIMULUS
for i = 1:nfilts
    temp = filt(i,:,:);
%     filt_stim(i,:) = zscore(stim_emb*temp(:));
%     xvfilt_stim(i,:) = zscore(xvstim_emb*temp(:));
    filt_stim(i,:) = stim_emb*temp(:);
    xvfilt_stim(i,:) = xvstim_emb*temp(:);
end


%% create spike function
%pass through internal NLs
beta = 3; %2.5
theta = 0;
coefs = ones(1,nfilts).*cvals;
g = zeros(1,NT);
xvg = zeros(1,xvNT);
for i = 1:nfilts
    lfilt_stim(i,:) = 1/beta*log(1+exp(beta*(filt_stim(i,:)-theta)));
    xvlfilt_stim(i,:) = 1/beta*log(1+exp(beta*(xvfilt_stim(i,:)-theta)));
    g = g + coefs(i)*lfilt_stim(i,:);
    xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
end
sg = std(g);
g = g/sg;
xvg = xvg/sg;

target_rate = 50; %in Hz
target_pspike = target_rate*dt;
cur_theta = 1.75; %1.5  2.5
cur_beta = 4; %4  4
p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
% p_spike = (g-cur_theta).^2; p_spike(p_spike<0) = 0;
% p_spike = exp((g-1.5)/1);
xvp_spike = 1/cur_beta*log(1+exp(cur_beta*(xvg-cur_theta)));
scale_f = target_rate/(mean(p_spike)/dt);
p_spike = p_spike*scale_f;
% p_spike(p_spike/dt > 500) = 500*dt;
xvp_spike = xvp_spike*scale_f;
% xvp_spike(xvp_spike/dt > 500) = 500*dt;


%%
close all
spikes = poissrnd(p_spike);
rbins = (find(spikes>0.5));
nsp = spikes(rbins);
spk_vals = unique(spikes); spk_vals(spk_vals==0) = [];
spikebins = [];
for i = 1:length(spk_vals)
    cur_set = find(spikes == spk_vals(i));
    spikebins = [spikebins; repmat(cur_set(:),spk_vals(i),1)];
end
% spikebins = unique(spikebins); %max 1 spike per bin
fprintf('Nspks: %d\n',length(spikebins));

xvspikes = poissrnd(xvp_spike);
xvrbins = (find(xvspikes>0.5));
xvnsp = spikes(xvrbins);
xvspk_vals = unique(xvspikes); xvspk_vals(xvspk_vals==0) = [];
xvspikebins = [];
for i = 1:length(xvspk_vals)
    cur_set = find(xvspikes == xvspk_vals(i));
    xvspikebins = [xvspikebins; repmat(cur_set(:),xvspk_vals(i),1)];
end


%% QUAD MOD
sdim = SDIM; flen = 14;
stim_params.spatial_dims = 1;
stim_params.sdim = sdim;
stim_params.flen = flen; 

lam_vals = 0;
defmod.lambda_d2XT = 0;
nmods = [4];
defmod.lambda_L1x = lam_vals;

init_kerns = randn(flen*sdim,nmods);
init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
init_signs = ones(1,nmods);
clear kern_types
kern_types{1} = 'lin';
for j = 2:nmods
    kern_types{j} = 'quad';
end
quad_mod = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
quad_mod = fitGNM_filters(quad_mod,stim_emb,spikebins,'none',[],1e-4,1e-6);
[~, ~, ~, ~, g] = getLL_GNM(quad_mod,stim_emb,spikebins,'none');
quad_mod = fitGNM_spkNL(quad_mod,g,spikebins,0);
quad_xvLL= getLL_GNM(quad_mod,xvstim_emb,xvspikebins,'none')

cur_xvLL = quad_xvLL;
dxvLL = Inf;
cur_mod = quad_mod;
while dxvLL > 1e-4
    quad_modr = fitGNM_filters(cur_mod,stim_emb,spikebins,'none',[],1e-4,1e-6);
    [~, ~, ~, ~, g] = getLL_GNM(quad_modr,stim_emb,spikebins,'none');
    quad_modr = fitGNM_spkNL(quad_modr,g,spikebins,0);
    cur_mod = quad_modr;
    quadr_xvLL= getLL_GNM(quad_modr,xvstim_emb,xvspikebins,'none')
    dxvLL = cur_xvLL-quadr_xvLL
    cur_xvLL = quadr_xvLL;
end

%%
sdim = SDIM; flen = 14;
stim_params.spatial_dims = 1;
stim_params.sdim = sdim;
stim_params.flen = flen;

defmod.lambda_d2XT = 0;

clear gnmr*
lam_vals = 0;
defmod.lambda_L1x = lam_vals;
nmods = 6;
init_kerns = randn(flen*sdim,nmods);

init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
init_signs = ones(1,nmods);

clear kern_types
for j = 1:nmods
    kern_types{j} = 'threshlin';
end
gnm = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
gnm = fitGNM_filters(gnm,stim_emb,spikebins,'none',[],1e-4,1e-6);
[~, ~, ~, ~, g] = getLL_GNM(gnm,stim_emb,spikebins,'none');
gnm = fitGNM_spkNL(gnm,g,spikebins,0);
gnm_xvLL = getLL_GNM(gnm,xvstim_emb,xvspikebins,'none')

cur_xvLL = gnm_xvLL;
dxvLL = Inf;
cur_mod = gnm;
while dxvLL > 1e-4
    gnmr = fitGNM_filters(cur_mod,stim_emb,spikebins,'none',[],1e-4,1e-6);
    [~, ~, ~, ~, g] = getLL_GNM(gnmr,stim_emb,spikebins,'none');
    gnmr = fitGNM_spkNL(gnmr,g,spikebins,0);
    cur_mod = gnmr;
    gnmr_xvLL= getLL_GNM(gnmr,xvstim_emb,xvspikebins,'none')
    dxvLL = cur_xvLL-gnmr_xvLL
    cur_xvLL = gnmr_xvLL;
end

[~,coms,peak_locs] = get_filter_coms_1d(gnm);
[~,ord] = sort(coms);
gnmr.mods = gnmr.mods(ord);

%for fitting internal NLs
gnmr = adjust_all_reg(gnmr,'lnl2',1000); %1000
gnmr = adjust_all_reg(gnmr,'nltype','uncon');
gnmr = adjust_all_reg(gnmr,'nlmon',1);
gnmr = fitGNM_internal_NLs(gnmr,stim_emb,spikebins,0,0);
[~, ~, ~, ~, g] = getLL_GNM(gnmr,stim_emb,spikebins,'none');
gnmr = fitGNM_spkNL(gnmr,g,spikebins,0);
gnmr_xvLL = getLL_GNM(gnmr,xvstim_emb,xvspikebins,'none')
gnmr2 = fitGNM_filters(gnmr,stim_emb,spikebins,'none',[],1e-4,1e-6);
gnmr2 = fitGNM_internal_NLs(gnmr2,stim_emb,spikebins,0,0);
[~, ~, ~, ~, g] = getLL_GNM(gnmr2,stim_emb,spikebins,'none');
gnmr2 = fitGNM_spkNL(gnmr2,g,spikebins,0);
gnmr2_xvLL = getLL_GNM(gnmr2,xvstim_emb,xvspikebins,'none')

