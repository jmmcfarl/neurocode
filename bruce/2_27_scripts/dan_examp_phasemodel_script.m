cd ~/Data/bruce/2_27_12/M232
load ./random_bar_phase_models_lfp_nsm2
load ./random_bar_phasemods_stimmods

%%
%SUS are [3 4 5 14 15 17]
%DUPLICATES: [22 23 24] [11 12] and maybe [2 3]

close all
for cc = 1:24
    fprintf('Stim only LL imp: %.4f bits/spk\n',(xv_null_LL(cc) - xv_so_LL(cc))/log(2));
    fprintf('Sac + Stim LL imp: %.4f bits/spk\n',(xv_null_LL(cc) - xv_sac_LL(cc))/log(2));
    fprintf('Stim + phase + sac LL imp: %.4f bits/spk\n',(xv_null_LL(cc) - xv_post_LL(cc))/log(2));
    
    %% PLOT SACCADE MODULATION KERNELS
    subplot(3,2,1);hold on
    plot(tent_centers*new_dt,so_sac_kern(cc,:),'b') %all saccades pooled
    hold on
    plot(tent_centers*new_dt,so_msac_kern(cc,:),'r') %microsaccades
    plot(tent_centers*new_dt,so_fsac_kern(cc,:),'k') %outgoing saccades
    plot(tent_centers*new_dt,so_ssac_kern(cc,:),'color',[0.2 0.8 0.2]) %returning saccades
    axis tight
    xlim([-0.2 0.4])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')
    
    %% PLOT STIMULUS STRF
    subplot(3,2,2);
    k_mat = glm_fit(cc).mods(1).k;
    imagesc(un_bar_pos,(0:-1:-flen+1)*new_dt,flipud(reshape(k_mat,flen,n_bar_pos)));
    cca = max(abs(k_mat(:)));
    caxis([-cca cca]*0.85);
    xlabel('Bar position deg)')
    ylabel('Time lag (s)')
    
    %% PLOT LFP MODEL AMPLITUDE KERNEL
    cur_ampkern = po_ampkern(cc,:); %without stimulus
    %     cur_ampkern = post_ampkern(cc,:); %with stimulus
    ca1 = max(cur_ampkern);
    cur_ampkern = [reshape(cur_ampkern',length(wfreqs),24)'; zeros(1,length(wfreqs))]; %pad with a zero row for pcolor plotting
    subplot(3,2,3)
    pcolor(wfreqs,1:25,cur_ampkern);
    shading flat; %doesn't look as good, but no artifacts
    shading interp; %notice interpolation artifact at the boundaries (specifically low-freq boundary)
    set(gca,'xscale','log');
    caxis([0 ca1]*0.75);
    xlabel('Frequency (Hz)')
    ylabel('Channel')
    
    %% PLOT PHASE KERNEL
    cur_phase_kern = po_phasekern(cc,:)'; %without stim
    %     cur_phase_kern = post_phasekern(cc,:)'; %with stim
    cur_phase_kern = [reshape(cur_phase_kern',length(wfreqs),24)'; zeros(1,length(wfreqs))]; %pad with a zero row for pcolor plotting
    subplot(3,2,4)
    pcolor(wfreqs,1:25,cur_phase_kern);
    shading flat;
    set(gca,'xscale','log');
    caxis([0 360]);
    title('Phase kernel')
    xlabel('Frequency (Hz)')
    ylabel('Channel')
    
    %%  PLOT SINE AND COSINE KERNELS
    cur_sine_kern = po_sinphase_sfilt(cc,:);
    cur_cos_kern = po_sinphase_cfilt(cc,:);
    %     cur_sine_kern = post_sinphase_sfilt(cc,:);
    %     cur_cos_kern = post_sinphase_cfilt(cc,:);
    
    ca = max([max(abs(cur_cos_kern)) max(abs(cur_sine_kern))]);
    
    cur_cos_kern = [reshape(cur_cos_kern',length(wfreqs),24)'; zeros(1,length(wfreqs))]; %pad with a zero row for pcolor plotting
    cur_sine_kern = [reshape(cur_sine_kern',length(wfreqs),24)'; zeros(1,length(wfreqs))]; %pad with a zero row for pcolor plotting
    
    subplot(3,2,5)
    pcolor(wfreqs,1:25,cur_sine_kern);
    shading interp;
    set(gca,'xscale','log');
    caxis(0.8*[-ca ca])
    title('Sine kernel')
    xlabel('Frequency (Hz)')
    ylabel('Depth (mm)')
    
    subplot(3,2,6)
    pcolor(wfreqs,1:25,cur_cos_kern);
    shading interp;
    set(gca,'xscale','log');
    caxis(0.8*[-ca ca])
    title('Cosine kernel')
    xlabel('Frequency (Hz)')
    ylabel('Depth (mm)')
    
    cc
    pause
    clf
end

