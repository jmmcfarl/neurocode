clear all
% close all
fig_dir = '/home/james/Desktop/';

dt = 0.01;

dx = 0.02;
xax = -0.4:dx:0.4;

tf = 4;
rf_mean = 0;
rf_psf = 5;
rf_std = 1/rf_psf/2.5;
gauss = exp(-(xax-rf_mean).^2/(2*rf_std^2));

n_rpts = 10;
NT = round(100/dt);

rel_quad = 0.;

delta_t = 0.01;

poss_ES = [0.1 0.2 0.4 0.6 0.8 1 1.5 2];
base_eye_SD = 0.1;
drift_fract = 0.075;
%%
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
%%
taxis = 0:dt:(NT-1)*dt;

sine_phase = 0:dt*tf*2*pi:(NT-1)*dt*tf*2*pi;
sine_stim = sin(bsxfun(@plus,xax*rf_psf*2*pi,sine_phase'));

NFFT = 2^nextpow2(NT);
f = 1/dt/2*linspace(0,1,NFFT/2+1);
[~,freqloc] = min(abs(f-tf));

%%
stim_params = NMMcreate_stim_params([flen length(xax)],dt);
X = create_time_embedding(sine_stim,stim_params);
%%
filt_out = X*RF(:);
qfilt_out = (X*RF(:)).^2 + (X*RFq(:)).^2;
tot_filt_out = filt_out*(1-rel_quad) + qfilt_out*rel_quad;
%%
spk_thresh = -1.5; %-1.5
spk_beta = 1; %1
spk_offset = 0;
mean_spk_rate = 50*dt;
out_rate = log(1+exp(spk_beta*tot_filt_out + spk_thresh)) + spk_offset;
sc_fac = mean_spk_rate/mean(out_rate);
out_rate = out_rate*sc_fac;

null_rate = (log(1+exp(spk_thresh))+spk_offset)*sc_fac;
%%
for nn = 1:n_rpts
    fprintf('Rep %d of %d\n',nn,n_rpts);
    
    out_cnts = poissrnd(out_rate); %sample spikes for perfect fixation condition
    
    P1 = fft(out_rate,NFFT)/NT;
    P_true(nn,:) = 2*abs(P1(1:NFFT/2+1)); %single-sided amplitude spectrum of rates
    
    C1 = fft(out_cnts,NFFT)/NT;
    C_true(nn,:) = 2*abs(C1(1:NFFT/2+1)); %s-sided amp spectrum of counts
    
    %F1/F0, note this uses the peak-to-peak amp of the sinusoide at best
    %freq
    simp_ind(nn) = 2*P_true(nn,freqloc)/mean(out_rate);
    simp_indC(nn) = 2*C_true(nn,freqloc)/mean(out_rate);
    
    for ee = 1:length(poss_ES)
        eye_std = poss_ES(ee)*base_eye_SD;
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
        
        %% CREATE SIM EYE POSITION CONSISTING ONLY OF MICROSACS (intantaneous saccades)
        sim_eyepos = nan(NT,1);
        for ii = 1:length(fix_start_inds)-1
            cur_eye_rnd = randn*eye_std;
            cur_inds = fix_start_inds(ii):(fix_start_inds(ii+1));
            sim_eyepos(cur_inds) = cur_eye_rnd;
        end
        cur_eye_rnd = randn*eye_std;
        cur_inds = fix_start_inds(end):NT;
        sim_eyepos(cur_inds) = cur_eye_rnd;
        
        %% construct retinal stimulus
        sim_eyepos_phase = sim_eyepos*rf_psf*2*pi;
        eye_sine_phase = sine_phase' + sim_eyepos_phase;
        eye_sine_stim = sin(bsxfun(@plus,xax*rf_psf*2*pi,eye_sine_phase));
        
        %% output of RF
        X = create_time_embedding(eye_sine_stim,stim_params);
        filt_out = X*RF(:);
 qfilt_out = (X*RF(:)).^2 + (X*RFq(:)).^2;
tot_filt_out = filt_out*(1-rel_quad) + qfilt_out*rel_quad;
       %% predicted rate
        eye_out_rate = log(1+exp(spk_beta*tot_filt_out + spk_thresh))+spk_offset;
        eye_out_rate = eye_out_rate*sc_fac;
        eye_out_cnts = poissrnd(eye_out_rate);
        
        %% now construct sim eye positoin with microsacs and drift (random-walk)
        rwalk_std = eye_std*drift_fract;
        sim_eyepos_drift = nan(NT,1);
        for ii = 1:length(fix_start_inds)-1
            cur_eye_rnd = randn*eye_std;
            cur_inds = fix_start_inds(ii):(fix_start_inds(ii+1));
            rwalk = cumsum(randn(length(cur_inds),1)*rwalk_std);
            sim_eyepos_drift(cur_inds) = cur_eye_rnd + rwalk;
        end
        cur_eye_rnd = randn*eye_std;
        cur_inds = fix_start_inds(end):NT;
        rwalk = cumsum(randn(length(cur_inds),1)*rwalk_std);
        sim_eyepos_drift(cur_inds) = cur_eye_rnd + rwalk;
        
        %% construct corresponding retinal stim
        sim_eyepos_phase_drift = sim_eyepos_drift*rf_psf*2*pi;
        eye_sine_phase_drift = sine_phase' + sim_eyepos_phase_drift;
        eye_sine_stim_drift = sin(bsxfun(@plus,xax*rf_psf*2*pi,eye_sine_phase_drift));
        
        %% output of RF
        X = create_time_embedding(eye_sine_stim_drift,stim_params);
        filt_out = X*RF(:);
qfilt_out = (X*RF(:)).^2 + (X*RFq(:)).^2;
tot_filt_out = filt_out*(1-rel_quad) + qfilt_out*rel_quad;
        %%
        eye_out_rate_drift = log(1+exp(spk_beta*tot_filt_out + spk_thresh))+spk_offset;
        eye_out_rate_drift = eye_out_rate_drift*sc_fac;
        eye_out_cnts_drift = poissrnd(eye_out_rate_drift);
        
        %% compute power spectra of output rate and counts
        Ptemp = fft(eye_out_rate,NFFT)/NT;
        Peye(nn,ee,:) = 2*abs(Ptemp(1:NFFT/2+1));
        Ptemp = fft(eye_out_rate_drift,NFFT)/NT;
        Peyed(nn,ee,:) = 2*abs(Ptemp(1:NFFT/2+1));
        
        Ptemp = fft(eye_out_cnts,NFFT)/NT;
        Ceye(nn,ee,:) = 2*abs(Ptemp(1:NFFT/2+1));
        Ptemp = fft(eye_out_cnts_drift,NFFT)/NT;
        Ceyed(nn,ee,:) = 2*abs(Ptemp(1:NFFT/2+1));
        %% compute simplicity and mean rates
        simp_ind_eye(nn,ee) = 2*Peye(nn,ee,freqloc)/mean(eye_out_rate);
        simp_ind_eyed(nn,ee) = 2*Peyed(nn,ee,freqloc)/mean(eye_out_rate_drift);
        
        simp_indC_eye(nn,ee) = 2*Ceye(nn,ee,freqloc)/mean(eye_out_rate);
        simp_indC_eyed(nn,ee) = 2*Ceyed(nn,ee,freqloc)/mean(eye_out_rate_drift);
        
        mrate_eye(nn,ee) = mean(eye_out_rate);
        mrate_eyed(nn,ee) = mean(eye_out_rate_drift);
    end
end
%%
corr_PSF = rf_psf.*poss_ES;

if n_rpts > 1
    h1 = figure(); hold on
    errorbar(corr_PSF,mean(simp_indC_eyed),std(simp_indC_eyed),'o-')
    errorbar(corr_PSF,mean(simp_indC_eye),std(simp_indC_eye),'ko-')
    line(corr_PSF([1 end]),[mean(simp_indC) mean(simp_indC)],'color','r');
    line(corr_PSF([1 end]),[mean(simp_indC) mean(simp_indC)]-[std(simp_indC) std(simp_indC)],'color','r','linestyle','--');
    line(corr_PSF([1 end]),[mean(simp_indC) mean(simp_indC)]+[std(simp_indC) std(simp_indC)],'color','r','linestyle','--');
    xlabel('Preferred SF (cyc/deg)');
    ylabel('F1/F0');
    figufy(h1);
else
    h1 = figure(); hold on
    plot(corr_PSF,simp_indC_eyed,'o-')
    line(corr_PSF([1 end]),[simp_indC simp_indC],'color','r');
    xlabel('Preferred SF (cyc/deg)');
    ylabel('F1/F0');
    figufy(h1);
end

%%
fig_width = 4;
rel_height = 0.8;
% fname = [fig_dir 'simp_EM.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);


%%
% figure; hold on
% errorbar(poss_ES*5,mean(mrate_eye),std(mrate_eye),'ko-')
% errorbar(poss_ES*5,mean(mrate_eyed),std(mrate_eyed),'o-')
% line(poss_ES([1 end])*5,[0.5 0.5],'color','r');

%%
% sm_fac = 0;
% % nn = 1;
% ee = 6
% h2 = figure();
% % subplot(2,1,1)
% % plot(f,smooth(C1,sm_fac),'k')
% % plot(f,C1,'k')
% shadedErrorBar(f,squeeze(mean(C_true)),squeeze(std(C_true)),'k')
% hold on
% % cmap = jet(length(poss_ES));
% % for ee = 1:length(poss_ES)
% % plot(f,smooth(squeeze(sqrt(Peye(nn,ee,:))),5),'r')
% % plot(f,smooth(squeeze(sqrt(Ceyed(nn,ee,:))),sm_fac),'r')
% % shadedErrorBar(f,squeeze(mean(smooth(sqrt(Ceyed(:,ee,:)),sm_fac))),squeeze(std(smooth(sqrt(Ceyed(:,ee,:)),sm_fac))),'r')
% % shadedErrorBar(f,squeeze(mean(smooth(Ceyed(:,ee,:),sm_fac))),squeeze(std(smooth(Ceyed(:,ee,:),sm_fac))),'r')
% shadedErrorBar(f,squeeze(mean(Ceyed(:,ee,:))),squeeze(std(Ceyed(:,ee,:))),'r')
% legend('No EM','With EM');
% % end
% xlim([0 15])
%
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% figufy(h2);


%%

fig_width = 4;
rel_height = 0.8;
% fname = [fig_dir 'simp_EM_pspec.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
% h3 = figure();
% imagesc(xax,dt/2:dt:(flen*dt-dt/2),RF);colormap(gray);set(gca,'ydir','normal')
% xlabel('Position (deg)');
% ylabel('Lag (s)');
%
% fig_width = 4;
% rel_height = 0.8;
% % fname = [fig_dir 'simp_RF.pdf'];
% % exportfig(h3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
% nn = 1;
% h4 = figure();
% plot(taxis,out_rate/dt,'k')
% hold on
% plot(taxis,eye_out_rate_drift/dt,'r')
% xlim([4 6])
% legend('No EM','EM');
% xlabel('Time (s)')
% ylabel('Firing rate (Hz)');

%%
% fig_width = 6;
% rel_height = 0.75;
% figufy(h4)
% fname = [fig_dir 'simp_EM_examprate.pdf'];
% exportfig(h4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
