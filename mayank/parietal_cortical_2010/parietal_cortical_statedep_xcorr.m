
%used_ups: used up transition indices
%used_downs: used down transition indices
%used_ups8
%used_downs8

clear all
close all

load F:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
load F:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8
addpath('F:\Code\Chronux\spectral_analysis\continuous\')
addpath('F:\WC_Germany\parietal_cortical_2010\')

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = 2016/2;
lcf1 = 4;
hcf1 = 10;
lcf2 = 10;
hcf2 = 20;
lcf3 = 20;
hcf3 = 40;
maxlag = round(Fsd*0.4);
lags = -maxlag:maxlag;

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

sess_data = sess_data(parietal);
desynch_start_times = desynch_start_times(parietal);
desynch_stop_times = desynch_stop_times(parietal);
%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];

n = length(sess_data);
xc1_up = zeros(n,length(lags));
xc2_up = zeros(n,length(lags));
xc3_up = zeros(n,length(lags));
xc1_down = zeros(n,length(lags));
xc2_down = zeros(n,length(lags));
xc3_down = zeros(n,length(lags));
xc1_8up = zeros(n,length(lags));
xc2_8up = zeros(n,length(lags));
xc3_8up = zeros(n,length(lags));
xc1_8down = zeros(n,length(lags));
xc2_8down = zeros(n,length(lags));
xc3_8down = zeros(n,length(lags));

for d = 1:n
    
    to_dir = sess_data(d).directory;
    to_dir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(to_dir);
    
    load used_data wcv_minus_spike lf8
    lf81 = get_hf_signal(lf8,raw_Fs,Fsd,[lcf1 hcf1]);
    lf82 = get_hf_signal(lf8,raw_Fs,Fsd,[lcf2 hcf2]);
    lf83 = get_hf_signal(lf8,raw_Fs,Fsd,[lcf3 hcf3]);
    wcv1 = get_hf_signal(wcv_minus_spike,raw_Fs,Fsd,[lcf1 hcf1]);
    wcv2 = get_hf_signal(wcv_minus_spike,raw_Fs,Fsd,[lcf2 hcf2]);
    wcv3 = get_hf_signal(wcv_minus_spike,raw_Fs,Fsd,[lcf3 hcf3]);
    
    %% get used state transition times
    load hsmm_state_seq_lf_pert15
    load hsmm_state_seq8_lf_pert15
    mp_state_seq =  hsmm_bbstate_seq;
    old_t = (1:length(hsmm_bbstate_seq))/Fs_bb;
    new_t = (1:length(lf81))/Fsd;
    mp_state_seq = round(interp1(old_t,mp_state_seq,new_t));
    lfp_state_seq = hsmm_bbstate_seq8;
    lfp_state_seq = round(interp1(old_t,lfp_state_seq,new_t));
    
    up_trans = find(mp_state_seq(1:end-1)==1 & mp_state_seq(2:end)==2);
    down_trans = find(mp_state_seq(1:end-1)==2 & mp_state_seq(2:end)==1);
    down_trans(down_trans < up_trans(1)) = [];
    up_trans(up_trans > down_trans(end)) = [];
    bad_states = [];
    for i = 1:length(desynch_start_times{d})
        if up_trans*Fsd > desynch_start_times{d}(i) & up_trans*Fsd < desynch_stop_times{d}(i)
            bad_states = [bad_states i];
        end
    end
    up_trans(bad_states) = [];
    down_trans(bad_states) = [];
    
    up_trans8 = find(lfp_state_seq(1:end-1)==1 & lfp_state_seq(2:end)==2);
    down_trans8 = find(lfp_state_seq(1:end-1)==2 & lfp_state_seq(2:end)==1);
    down_trans8(down_trans8 < up_trans8(1)) = [];
    up_trans8(up_trans8 > down_trans8(end)) = [];
    bad_states = [];
    for i = 1:length(desynch_start_times{d})
        if up_trans8*Fsd > desynch_start_times{d}(i) & up_trans8*Fsd < desynch_stop_times{d}(i)
            bad_states = [bad_states i];
        end
    end
    up_trans8(bad_states) = [];
    down_trans8(bad_states) = [];
    
    up_starts = up_trans;
    up_stops = down_trans;
    down_starts = down_trans(1:end-1);
    down_stops = up_trans(2:end);
    up_starts8 = up_trans8;
    up_stops8 = down_trans8;
    down_starts8 = down_trans8(1:end-1);
    down_stops8 = up_trans8(2:end);
    
    t_per = round(0.15*Fsd);
    up_starts = up_starts + t_per;
    up_stops = up_stops - t_per;
    down_starts = down_starts + t_per;
    down_stops = down_stops - t_per;
    up_starts8 = up_starts8 + t_per;
    up_stops8 = up_stops8 - t_per;
    down_starts8 = down_starts8 + t_per;
    down_stops8 = down_stops8 - t_per;
    
    up_durs = (up_stops-up_starts)/Fsd;
    too_short = find(up_durs < 0.4);
    up_starts(too_short) = [];
    up_stops(too_short) = [];
    down_durs = (down_stops-down_starts)/Fsd;
    too_short = find(down_durs < 0.4);
    down_starts(too_short) = [];
    down_stops(too_short) = [];
    
    up_durs8 = (up_stops8-up_starts8)/Fsd;
    too_short = find(up_durs8 < 0.4);
    up_starts8(too_short) = [];
    up_stops8(too_short) = [];
    down_durs8 = (down_stops8-down_starts8)/Fsd;
    too_short = find(down_durs8 < 0.4);
    down_starts8(too_short) = [];
    down_stops8(too_short) = [];
    
    up_markers = [up_starts' up_stops'];
    down_markers = [down_starts' down_stops'];
    up_markers8 = [up_starts8' up_stops8'];
    down_markers8 = [down_starts8' down_stops8'];
    
    n_ups = length(up_starts);
    for i = 1:n_ups
        cur_seg = up_starts(i):up_stops(i);
        xc1_up(d,:) = xc1_up(d,:) + xcov(wcv1(cur_seg),lf81(cur_seg),maxlag,'coeff')';
        xc2_up(d,:) = xc2_up(d,:) + xcov(wcv2(cur_seg),lf82(cur_seg),maxlag,'coeff')';
        xc3_up(d,:) = xc3_up(d,:) + xcov(wcv3(cur_seg),lf83(cur_seg),maxlag,'coeff')';
    end
    xc1_up(d,:) = xc1_up(d,:)/n_ups;
    xc2_up(d,:) = xc2_up(d,:)/n_ups;
    xc3_up(d,:) = xc3_up(d,:)/n_ups;
    
    n_downs = length(down_starts);
    for i = 1:n_downs
        cur_seg = down_starts(i):down_stops(i);
        xc1_down(d,:) = xc1_down(d,:)  + xcov(wcv1(cur_seg),lf81(cur_seg),maxlag,'coeff')';
        xc2_down(d,:) = xc2_down(d,:)  + xcov(wcv2(cur_seg),lf82(cur_seg),maxlag,'coeff')';
        xc3_down(d,:) = xc3_down(d,:)  + xcov(wcv3(cur_seg),lf83(cur_seg),maxlag,'coeff')';
    end
    xc1_down(d,:) = xc1_down(d,:)/n_downs;
    xc2_down(d,:) = xc2_down(d,:)/n_downs;
    xc3_down(d,:) = xc3_down(d,:)/n_downs;
    
    n_8ups = length(up_starts8);
    for i = 1:n_8ups
        cur_seg = up_starts8(i):up_stops8(i);
        xc1_8up(d,:) = xc1_8up(d,:) + xcov(wcv1(cur_seg),lf81(cur_seg),maxlag,'coeff')';
        xc2_8up(d,:) = xc2_8up(d,:) + xcov(wcv2(cur_seg),lf82(cur_seg),maxlag,'coeff')';
        xc3_8up(d,:) = xc3_8up(d,:) + xcov(wcv3(cur_seg),lf83(cur_seg),maxlag,'coeff')';
    end
    xc1_8up(d,:) = xc1_8up(d,:)/n_ups;
    xc2_8up(d,:) = xc2_8up(d,:)/n_ups;
    xc3_8up(d,:) = xc3_8up(d,:)/n_ups;
    
    n_8downs = length(down_starts8);
    for i = 1:n_8downs
        cur_seg = down_starts8(i):down_stops8(i);
        xc1_8down(d,:) = xc1_8down(d,:)  + xcov(wcv1(cur_seg),lf81(cur_seg),maxlag,'coeff')';
        xc2_8down(d,:) = xc2_8down(d,:)  + xcov(wcv2(cur_seg),lf82(cur_seg),maxlag,'coeff')';
        xc3_8down(d,:) = xc3_8down(d,:)  + xcov(wcv3(cur_seg),lf83(cur_seg),maxlag,'coeff')';
    end
    xc1_8down(d,:) = xc1_8down(d,:)/n_downs;
    xc2_8down(d,:) = xc2_8down(d,:)/n_downs;
    xc3_8down(d,:) = xc3_8down(d,:)/n_downs;
    
end

cd F:\WC_Germany\parietal_cortical_2010\
save state_dep_xc_notrans xc* lags Fsd

%%
figure
subplot(2,1,1)
h=errorbar(lags/Fsd,mean(xc1_up),std(xc1_up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc1_down),std(xc1_down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
legend('MP up state','MP down state')
xlim([-maxlag maxlag]/Fsd)
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('4-10Hz','fontsize',18)
ylim([-0.15 0.15])
line([0 0],[-0.15 0.15],'color','k')
subplot(2,1,2)
h=errorbar(lags/Fsd,mean(xc1_8up),std(xc1_8up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc1_8down),std(xc1_8down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
legend('LFP up state','LFP down state')
xlim([-maxlag maxlag]/Fsd)
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('4-10Hz','fontsize',18)
ylim([-0.15 0.15])
line([0 0],[-0.15 0.15],'color','k')

figure
subplot(2,1,1)
h=errorbar(lags/Fsd,mean(xc2_up),std(xc2_up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc2_down),std(xc2_down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
legend('MP up state','MP down state')
xlim([-0.2 0.2])
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('10-20Hz','fontsize',18)
ylim([-0.1 0.1])
line([0 0],[-0.1 0.1],'color','k')
subplot(2,1,2)
h=errorbar(lags/Fsd,mean(xc2_8up),std(xc2_8up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc2_8down),std(xc2_8down)/sqrt(n),'r');
legend('LFP up state','LFP down state')
errorbar_tick(h,.001,'units')
xlim([-0.2 0.2])
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('10-20Hz','fontsize',18)
ylim([-0.1 0.1])
line([0 0],[-0.1 0.1],'color','k')

figure
subplot(2,1,1)
h=errorbar(lags/Fsd,mean(xc3_up),std(xc3_up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc3_down),std(xc3_down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
xlim([-0.1 0.1])
legend('MP up state','MP down state')
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('20-40Hz','fontsize',18)
ylim([-0.06 0.06])
line([0 0],[-0.06 0.06],'color','k')
subplot(2,1,2)
h=errorbar(lags/Fsd,mean(xc3_8up),std(xc3_8up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc3_8down),std(xc3_8down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
legend('LFP up state','LFP down state')
xlim([-0.1 0.1])
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('20-40Hz','fontsize',18)
ylim([-0.06 0.06])
line([0 0],[-0.06 0.06],'color','k')

%%
[pks_up1,pklcs_up1] = max(xc1_up,[],2);
[pks_up2,pklcs_up2] = max(xc2_up,[],2);
[pks_up3,pklcs_up3] = max(xc3_up,[],2);
[pks_down1,pklcs_down1] = max(xc1_down,[],2);
[pks_down2,pklcs_down2] = max(xc2_down,[],2);
[pks_down3,pklcs_down3] = max(xc3_down,[],2);

[pks_8up1,pklcs_8up1] = max(xc1_8up,[],2);
[pks_8up2,pklcs_8up2] = max(xc2_8up,[],2);
[pks_8up3,pklcs_8up3] = max(xc3_8up,[],2);
[pks_8down1,pklcs_8down1] = max(xc1_8down,[],2);
[pks_8down2,pklcs_8down2] = max(xc2_8down,[],2);
[pks_8down3,pklcs_8down3] = max(xc3_8down,[],2);

figure
plot(lags(pklcs_up1)/Fsd,pks_up1,'.')
hold on
plot(lags(pklcs_down1)/Fsd,pks_down1,'r.')
xlabel('peak lag (s)','Fontsize',16)
ylabel('peak correlation','fontsize',16)

figure
plot(lags(pklcs_up2)/Fsd,pks_up2,'.')
hold on
plot(lags(pklcs_down2)/Fsd,pks_down2,'r.')
xlabel('peak lag (s)','Fontsize',16)
ylabel('peak correlation','fontsize',16)

figure
plot(lags(pklcs_up3)/Fsd,pks_up3,'.')
hold on
plot(lags(pklcs_down3)/Fsd,pks_down3,'r.')
xlabel('peak lag (s)','Fontsize',16)
ylabel('peak correlation','fontsize',16)

figure
plot(lags(pklcs_8up1)/Fsd,pks_8up1,'.')
hold on
plot(lags(pklcs_8down1)/Fsd,pks_8down1,'r.')
xlabel('peak lag (s)','Fontsize',16)
ylabel('peak correlation','fontsize',16)

figure
plot(lags(pklcs_8up2)/Fsd,pks_8up2,'.')
hold on
plot(lags(pklcs_8down2)/Fsd,pks_8down2,'r.')
xlabel('peak lag (s)','Fontsize',16)
ylabel('peak correlation','fontsize',16)

figure
plot(lags(pklcs_8up3)/Fsd,pks_8up3,'.')
hold on
plot(lags(pklcs_8down3)/Fsd,pks_8down3,'r.')
xlabel('peak lag (s)','Fontsize',16)
ylabel('peak correlation','fontsize',16)



%%
%%
figure
subplot(2,1,1)
h=errorbar(lags/Fsd,mean(xc1_up),std(xc1_up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc1_down),std(xc1_down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
legend('MP up state','MP down state')
xlim([-maxlag maxlag]/Fsd)
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('4-10Hz','fontsize',18)
ylim([-0.15 0.15])
line([0 0],[-0.15 0.15],'color','k')
subplot(2,1,2)
h=errorbar(lags/Fsd,mean(xc1_8up),std(xc1_8up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc1_8down),std(xc1_8down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
legend('LFP up state','LFP down state')
xlim([-maxlag maxlag]/Fsd)
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('4-10Hz','fontsize',18)
ylim([-0.15 0.15])
line([0 0],[-0.15 0.15],'color','k')

figure
subplot(2,1,1)
h=errorbar(lags/Fsd,mean(xc2_up),std(xc2_up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc2_down),std(xc2_down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
legend('MP up state','MP down state')
xlim([-0.2 0.2])
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('10-20Hz','fontsize',18)
ylim([-0.1 0.1])
line([0 0],[-0.1 0.1],'color','k')
subplot(2,1,2)
h=errorbar(lags/Fsd,mean(xc2_8up),std(xc2_8up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc2_8down),std(xc2_8down)/sqrt(n),'r');
legend('LFP up state','LFP down state')
errorbar_tick(h,.001,'units')
xlim([-0.2 0.2])
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('10-20Hz','fontsize',18)
ylim([-0.1 0.1])
line([0 0],[-0.1 0.1],'color','k')

figure
subplot(2,1,1)
h=errorbar(lags/Fsd,mean(xc3_up),std(xc3_up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc3_down),std(xc3_down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
xlim([-0.1 0.1])
legend('MP up state','MP down state')
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('20-40Hz','fontsize',18)
ylim([-0.06 0.06])
line([0 0],[-0.06 0.06],'color','k')
subplot(2,1,2)
h=errorbar(lags/Fsd,mean(xc3_8up),std(xc3_8up)/sqrt(n)); hold on
errorbar_tick(h,.001,'units')
h=errorbar(lags/Fsd,mean(xc3_8down),std(xc3_8down)/sqrt(n),'r');
errorbar_tick(h,.001,'units')
legend('LFP up state','LFP down state')
xlim([-0.1 0.1])
xlabel('Time (s)','fontsize',16)
ylabel('correlation','fontsize',16)
title('20-40Hz','fontsize',18)
ylim([-0.06 0.06])
line([0 0],[-0.06 0.06],'color','k')

%%
