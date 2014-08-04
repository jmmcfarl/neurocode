clear all
close all
fig_dir = '/home/james/Desktop/';

dt = 0.01;

dx = 0.02;
xax = -0.4:dx:0.4;

rf_mean = 0;
rf_psf = 5;
rf_std = 1/rf_psf/2.5;
gauss = exp(-(xax-rf_mean).^2/(2*rf_std^2));
tf = 4;
n_rpts = 2;
poss_tf = [0.5 1 2 4 8 12];
NT = round(300/dt);

% base_eye_SD = 0.1;
% drift_fract = 0.075;
base_eye_SD = 0.1;
drift_fract = 0.075;
%%
flen = 14;
RF = nan(flen,length(xax));
% rf_phase = 0:dt*tf*2*pi:flen*tf*2*pi;
rf_phase = 0:-dt*tf*2*pi:-flen*tf*2*pi;
amp_vec = gampdf((1:flen-1)-1,3,1.3);
amp_vec = [0 amp_vec];
for ii = 1:flen
    RF(ii,:) = amp_vec(ii)*gauss.*sin(2*pi*xax*rf_psf + rf_phase(ii));
end
%%
taxis = 0:dt:(NT-1)*dt;


%%
params.tapers = [5 9];
params.Fs = 1/dt;
%%
for nn = 1:n_rpts
    for ee = 1:length(poss_tf)
        tf = poss_tf(ee);
% rf_psf = poss_tf(ee);
        ee
        
        sine_phase = 0:dt*tf*2*pi:(NT-1)*dt*tf*2*pi;
        sine_stim = sin(bsxfun(@plus,xax*rf_psf*2*pi,sine_phase'));
        
        %%
        stim_params = NMMcreate_stim_params([flen length(xax)],dt);
        X = create_time_embedding(sine_stim,stim_params);
        %%
        filt_out = X*RF(:);
        %%
        spk_thresh = -1.5; %-1.5
        spk_beta = 1; %1
        spk_offset = 0;
        mean_spk_rate = 50*dt;
        out_rate = log(1+exp(spk_beta*filt_out + spk_thresh)) + spk_offset;
        sc_fac = 1;
        out_rate = out_rate*sc_fac;
        out_cnts = poissrnd(out_rate);
        
        [P1,f] = mtspectrumc(out_rate,params);
        [C1,f] = mtspectrumc(out_cnts,params);
        [~,freqloc] = min(abs(f-tf));
        simp_ind = sqrt(P1(freqloc))/mean(out_rate);
        simp_indC = sqrt(C1(freqloc))/mean(out_rate);
        
        
        %%
        eye_std = 0.1;
        fix_alpha = 2;
        fix_beta = 0.15;
        min_fix_dur = 0.05;
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
        sim_eyepos(isnan(sim_eyepos)) = randn*eye_std;
        
        %%
        sim_eyepos_phase = sim_eyepos*rf_psf*2*pi;
        eye_sine_phase = sine_phase' + sim_eyepos_phase;
        eye_sine_stim = sin(bsxfun(@plus,xax*rf_psf*2*pi,eye_sine_phase));
        
        %%
        X = create_time_embedding(eye_sine_stim,stim_params);
        %%
        filt_out = X*RF(:);
        %%
        eye_out_rate = log(1+exp(spk_beta*filt_out + spk_thresh))+spk_offset;
        eye_out_rate = eye_out_rate*sc_fac;
        eye_out_cnts = poissrnd(eye_out_rate);
        
        %%
        rwalk_std = eye_std*drift_fract;
        sim_eyepos_drift = nan(NT,1);
        for ii = 1:length(fix_start_inds)-1
            cur_eye_rnd = randn*eye_std;
            cur_inds = fix_start_inds(ii):(fix_start_inds(ii+1));
            rwalk = cumsum(randn(length(cur_inds),1)*rwalk_std);
            sim_eyepos_drift(cur_inds) = cur_eye_rnd + rwalk;
        end
        sim_eyepos_drift(isnan(sim_eyepos_drift)) = randn*eye_std;
        
        %%
        sim_eyepos_phase_drift = sim_eyepos_drift*rf_psf*2*pi;
        eye_sine_phase_drift = sine_phase' + sim_eyepos_phase_drift;
        eye_sine_stim_drift = sin(bsxfun(@plus,xax*rf_psf*2*pi,eye_sine_phase_drift));
        
        %%
        X = create_time_embedding(eye_sine_stim_drift,stim_params);
        %%
        filt_out = X*RF(:);
        %%
        eye_out_rate_drift = log(1+exp(spk_beta*filt_out + spk_thresh))+spk_offset;
        eye_out_rate_drift = eye_out_rate_drift*sc_fac;
        eye_out_cnts_drift = poissrnd(eye_out_rate_drift);
        
        %%
        [Peye(nn,ee,:),f] = mtspectrumc(eye_out_rate,params);
        [Peyed(nn,ee,:),f] = mtspectrumc(eye_out_rate_drift,params);
        
        [Ceye(nn,ee,:),f] = mtspectrumc(eye_out_cnts,params);
        [Ceyed(nn,ee,:),f] = mtspectrumc(eye_out_cnts_drift,params);
        
        %%
        simp_ind_eye(nn,ee) = sqrt(Peye(nn,ee,freqloc))/mean(eye_out_rate);
        simp_ind_eyed(nn,ee) = sqrt(Peyed(nn,ee,freqloc))/mean(eye_out_rate_drift);
        
        simp_indC_eye(nn,ee) = sqrt(Ceye(nn,ee,freqloc))/mean(eye_out_rate);
        simp_indC_eyed(nn,ee) = sqrt(Ceyed(nn,ee,freqloc))/mean(eye_out_rate_drift);
        
        mrate_eye(nn,ee) = mean(eye_out_rate);
        mrate_eyed(nn,ee) = mean(eye_out_rate_drift);
        mrate(nn,ee) = mean(out_rate);
    end
end
%%


%%
figure; hold on
errorbar(poss_tf,mean(mrate_eye),std(mrate_eye),'ko-')
errorbar(poss_tf,mean(mrate_eyed),std(mrate_eyed),'o-')
errorbar(poss_tf,mean(mrate),std(mrate),'ro-');

%%
sm_fac = 10;
% nn = 1;
ee = 6;
h2 = figure();
% subplot(2,1,1)
plot(f,smooth(sqrt(C1),sm_fac),'k')
hold on
% cmap = jet(length(poss_ES));
% for ee = 1:length(poss_ES)
% plot(f,smooth(squeeze(sqrt(Peye(nn,ee,:))),5),'r')
% plot(f,smooth(squeeze(sqrt(Ceyed(nn,ee,:))),sm_fac),'r')
shadedErrorBar(f,squeeze(mean(smooth(sqrt(Ceyed(:,ee,:)),sm_fac))),squeeze(std(smooth(sqrt(Ceyed(:,ee,:)),sm_fac))),'r')
legend('No EM','With EM');
% end
xlim([0 15])

xlabel('Frequency (Hz)');
ylabel('Amplitude');
figufy(h2);


%%

fig_width = 4;
rel_height = 0.8;
% fname = [fig_dir 'simp_EM_pspec.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
h3 = figure();
imagesc(xax,dt/2:dt:(flen*dt-dt/2),RF);colormap(gray);set(gca,'ydir','normal')
xlabel('Position (deg)');
ylabel('Lag (s)');

fig_width = 4;
rel_height = 0.8;
% fname = [fig_dir 'simp_RF.pdf'];
% exportfig(h3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
nn = 1;
h4 = figure();
plot(taxis,out_rate/dt,'k')
hold on
plot(taxis,eye_out_rate_drift/dt,'r')
xlim([4 6])
legend('No EM','EM');
xlabel('Time (s)')
ylabel('Firing rate (Hz)');

%%
fig_width = 6;
rel_height = 0.75;
figufy(h4)
% fname = [fig_dir 'simp_EM_examprate.pdf'];
% exportfig(h4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
