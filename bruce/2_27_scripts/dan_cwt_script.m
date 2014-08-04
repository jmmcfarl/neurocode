load ./examp_data
%% SPECIFY DESIRED WAVELET PSEUDOFREQUENCIES
nwfreqs = 50;
scales = logspace(log10(3),log10(100),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

%% COMPUTE WAVELET TRANSFORM
all_cwt_lfp = zeros(24,length(lfp_timed),length(wfreqs));
for cc = 1:24
    fprintf('Channel %d of %d\n',cc,24);
    temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
    all_cwt_lfp(cc,:,:) = temp';
end

%%
all_filt_bank = real(all_cwt_lfp); %treating the wavelet transform as a series of filter banks with center frequencies given by wfreqs
all_cwt_amp = abs(all_cwt_lfp); %magnitude of the wavelet transform (power given by squared amplitude, but amplitude is better for visualization)
all_cwt_phase =  angle(all_cwt_lfp); %phase of LFP at each pseudo-frequency

%%
t_range = [80 83];

target_eye_inds = find(eyets >= t_range(1) & eyets < t_range(2));
target_lfp_inds = find(lfp_timed >= t_range(1) & lfp_timed < t_range(2));

use_lfp_ch = 2;

f1 = figure;
subplot(2,1,1)
plot(eyets(target_eye_inds),eye_speed(target_eye_inds))
xlabel('Time (s)','fontsize',16)
ylabel('Eye speed (deg/s)','fontsize',16)
title('Eye speed','fontsize',18)
subplot(2,1,2)
plot(lfp_timed(target_lfp_inds),lfp_sampsd(target_lfp_inds,use_lfp_ch))
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude','fontsize',16)
title('LFP signal','fontsize',18)

f2 = figure;
subplot(2,1,1)
pcolor(lfp_timed(target_lfp_inds),wfreqs,squeeze(all_filt_bank(use_lfp_ch,target_lfp_inds,:))');shading interp
set(gca,'yscale','log');
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Wavelet transformed LFP (real component)','fontsize',18)
subplot(2,1,2)
pcolor(lfp_timed(target_lfp_inds),wfreqs,squeeze(all_cwt_amp(use_lfp_ch,target_lfp_inds,:))');shading interp
set(gca,'yscale','log')
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Wavelet amplitude','fontsize',18)
