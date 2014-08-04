clear all
% close all
% fig_dir = '/home/james/Analysis/bruce/ET_final/';
fig_dir = '/Volumes/james/Analysis/bruce/ET_final/';

% bar_width = 0.0283;
bar_width = 0.0565;
usfac = 4;
% dx = bar_width/4;
dt = 0.01;
flen = 12;

% xax = -1:dx:1;
pix_ax = -0.5:bar_width:0.5;
xax = linspace(pix_ax(1),pix_ax(end),length(pix_ax)*usfac);
dx = median(diff(xax));

NT = round(1000/dt);
taxis = 0:dt:(NT-1)*dt;
interp_taxis = 0:0.001:(NT-1)*dt;

rel_quad = 0;

rf_mean = 0;
rf_psf = 4;
rf_std = 1/rf_psf/2.2;
rf_phase = pi/2;
rf_amp = 1;
tf = 4;

eye_std = 0.0875;

gauss = exp(-(xax-rf_mean).^2/(2*rf_std^2));
flen = 14;
RF = nan(flen,length(xax));
RFq = nan(flen,length(xax));
rf_phase = 0:-dt*tf*2*pi:-flen*tf*2*pi;
amp_vec = gampdf((1:flen-1)-1,3,1.3);
amp_vec = [0 amp_vec];
for ii = 1:flen
    RF(ii,:) = amp_vec(ii)*gauss.*sin(2*pi*xax*rf_psf + rf_phase(ii));
    RFq(ii,:) = amp_vec(ii)*gauss.*sin(2*pi*xax*rf_psf + rf_phase(ii) + pi/2);
end

noise_frac = 0.1;
spk_noise_frac = 1.5;
%% GENERATE STIM
dds = 0.12;
Npix = length(pix_ax);
Xs = zeros(NT,Npix);
nzero = rand(NT,Npix) < dds;
polarity = rand(NT,Npix) > 0.5;
Xs(polarity & nzero) = 1;
Xs(nzero & ~polarity) = -1;

Xs = reshape(repmat(Xs,1,usfac),[NT Npix usfac]);
Xs = permute(Xs,[1 3 2]);
Xs = reshape(Xs,NT,length(xax));

stim_params = NMMcreate_stim_params([flen Npix*usfac],dt);
X = create_time_embedding(Xs,stim_params);
%%
spk_thresh = -2; %-2
spk_beta = 5;%5
mean_spk_rate = 20*dt;

filt_out = X*RF(:);
filtqp_out = (X*RFq(:)).^2 + (X*RF(:)).^2;
filt_out = (1-rel_quad)*filt_out + rel_quad*filtqp_out;

spk_noise = randn(size(filt_out))*std(filt_out)*spk_noise_frac;
filt_out = filt_out + spk_noise;

pred_rate = log(1+exp(spk_beta*filt_out + spk_thresh));
sc_fac = mean_spk_rate/mean(pred_rate);
pred_rate = pred_rate*sc_fac;
Robs = poissrnd(pred_rate);

sta = mean(bsxfun(@times,X,Robs));

%%
% poss_sc_facs = [0.25 0.5 1 2 4];
% poss_sc_facs = [0 0.1 0.2 0.4 0.6 0.8 1 1.5 2];
poss_sc_facs = [0 0.1 0.2 0.5 1 1.5 2];
for ee = 1:length(poss_sc_facs)
    ee
    eye_std = 0.1*poss_sc_facs(ee);
    fix_alpha = 2; %shape param
    fix_beta = 0.1875; %set to give mean fixation duration of about 375 ms (alpha*beta = mean)
    min_fix_dur = 0.05; %threshold minimum fixation duration
    fix_durs = gamrnd(fix_alpha,fix_beta,[1e5 1]);
    fix_durs(fix_durs < min_fix_dur) = [];
    
    fix_start_inds = round([1; cumsum(fix_durs/dt)]);
    extra = find(fix_start_inds > NT);
    if isempty(extra)
        error('Need more fixations')
    end
    fix_start_inds(extra) = [];
    
    sim_eyepos = nan(NT,1);
    for ii = 1:length(fix_start_inds)-1
        cur_eye_rnd = randn*eye_std;
        cur_inds = fix_start_inds(ii):(fix_start_inds(ii+1));
        sim_eyepos(cur_inds) = cur_eye_rnd;
    end
    cur_eye_rnd = randn*eye_std;
    cur_inds = fix_start_inds(end):NT;
    sim_eyepos(cur_inds) = cur_eye_rnd;
    
    max_shift = Npix*usfac*dx;
    sim_eyepos(sim_eyepos > max_shift) = max_shift;
    sim_eyepos(sim_eyepos < -max_shift) = -max_shift;
    %%
    shiftX = Xs;
    for ii = 1:NT
        shiftX(ii,:) = shift_matrix_Nd(shiftX(ii,:),-round(sim_eyepos(ii)/dx),2);
    end
    shiftX = create_time_embedding(shiftX,stim_params);
    sta_shifted = mean(bsxfun(@times,shiftX,Robs));
    
    %%
    noise = randn(size(sim_eyepos))*sqrt((var(sim_eyepos)*noise_frac));
    sim_eyepos = noise;
    shiftX_cor = Xs;
    for ii = 1:NT
        shiftX_cor(ii,:) = shift_matrix_Nd(shiftX_cor(ii,:),-round(sim_eyepos(ii)/dx),2);
    end
    shiftX_cor = create_time_embedding(shiftX_cor,stim_params);
    sta_shifted_cor = mean(bsxfun(@times,shiftX_cor,Robs));
    
    %%
    reg_params = NMMcreate_reg_params('lambda_d2XT',100);
    if rel_quad == 0
        mod1 = NMMinitialize_model(stim_params,1,{'lin'},reg_params);
    else
        mod1 = NMMinitialize_model(stim_params,[1 1 1],{'lin','quad','quad'},reg_params);
    end
    mod2 = mod1;
    
    mod1 = NMMfit_filters(mod1,Robs,shiftX_cor);
    mod2 = NMMfit_filters(mod2,Robs,shiftX);
    linfilt = mod1.mods(1).filtK;
    linfilt_shifted = mod2.mods(1).filtK;
    
    %%
    [~,~,outrate] = NMMmodel_eval(mod1,Robs,X);
    [~,~,outrate_rev] = NMMmodel_eval(mod1,Robs,-X);
    cv(ee) = mean(abs(outrate-outrate_rev))/mean(outrate);
    
    [~,~,outrate] = NMMmodel_eval(mod2,Robs,shiftX);
    [~,~,outrate_rev] = NMMmodel_eval(mod2,Robs,-shiftX);
    cv_before(ee) = mean(abs(outrate-outrate_rev))/mean(outrate);
    
end

%%
figure
plot(poss_sc_facs*rf_psf,cv_before,'o-')
hold on
plot(poss_sc_facs*rf_psf,cv,'ro-')