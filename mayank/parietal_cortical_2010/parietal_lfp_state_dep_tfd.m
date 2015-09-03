clear all
% close all

load G:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8_2_24_09
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

raw_Fs = 2016;
dsf = 16;
Fsd = raw_Fs/dsf;
params.Fs = Fsd;
params.tapers = [2 3];
params.pad = 1;
params.err = [2 0.05];
params.fpass = [3 40];
params.trialave = 1;
% movingwin = [0.2 0.05];
movingwin = [0.4 0.05];
niqf = 2016/2;
lcf = 3/niqf;
hcf = 45/niqf;
[b,a] = butter(2,[lcf hcf]);
back_wins = 30;
state_edges = [0 round(Fsd*0.15)];
maxlag = round(0.3*Fsd);

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times(interneurons) = [];


% %get rid of sessions without lfp
% bad_sessions = [];
% for i = 1:length(sess_data)
%     if ismember(0,sess_data(i).gains)
%         bad_sessions = [bad_sessions i];
%     end
% end
% sess_data(bad_sessions) = [];
% desynch_start_times(bad_sessions) = [];
% desynch_stop_times(bad_sessions) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);
% desynch_times = desynch_times(thomas_el);

n = length(sess_data);

for d = 1:n
    
    to_dir = sess_data(d).directory;
    to_dir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(to_dir);
    
    load used_data wcv_minus_spike lf8 lf5 lf4
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = downsample(wcv_f,dsf)/sess_data(d).gains(1);
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
%     lf5_f = filtfilt(b,a,lf5);
%     lf5_f = downsample(lf5_f,dsf)/sess_data(d).gains(5);
    if sess_data(d).thom_elec
        lf4_f = filtfilt(b,a,lf4);
        lf4_f = downsample(lf4_f,dsf)/sess_data(d).gains(4);
    end
        
    %% get used state transition times
%     load hsmm_state_seq_seg_lf_pert15
    load hsmm_state_seq8_seg_lf_pert15
%     mp_state_seq =  hsmm_bbstate_seq;
    lfp_state_seq = hsmm_bbstate_seq8;
        
    rec_edge_buffer = 2*Fsd;
    [up_segs,down_segs,n_upsegs{d},n_downsegs{d}] = get_parsed_state_segments_cohgram_seg...
        (lfp_state_seq,252,Fsd,movingwin,hmm8.UDS_segs,rec_edge_buffer,state_edges,back_wins);
    
    used_segs = down_segs;
    n_used_segs{d} = n_downsegs{d};
    %
    %require at least 10 states lasting this long for each cell
    used_segs(n_used_segs{d} < 10) = [];
    n_used_segs{d}(n_used_segs{d} < 10) = [];
    for n = 1:length(used_segs)
        up_markers1 = [used_segs{n}(:) used_segs{n}(:)+round(Fsd*movingwin(1))-1];
        up_markers2 = up_markers1;

%         %%% for random pairings
%         cn = length(used_segs{n});
%         up_markers2 = up_markers1(randperm(cn));
%         up_markers2 = [up_markers2(:) up_markers2(:) + round(Fsd*movingwin(1))-1];        
%         %%%
               
        [Cw8,Phimn,Smn,Sw8,f] = ...
            coherencyc_unequal_length_trials_jmm_fixedsegs([wcv_f lf8_f], movingwin, params, up_markers1, up_markers2);
%         [Cw5,Phimn,Smn,Sw5,f] = ...
%             coherencyc_unequal_length_trials_jmm_fixedsegs([wcv_f lf5_f], movingwin, params, up_markers1, up_markers2);
        if sess_data(d).thom_elec
            [Cw4,Phimn,Smn,Sw4,f] = ...
                coherencyc_unequal_length_trials_jmm_fixedsegs([wcv_f lf4_f], movingwin, params, up_markers1, up_markers2);
            [C84,Phimn,Smn,S84,f] = ...
                coherencyc_unequal_length_trials_jmm_fixedsegs([lf8_f lf4_f], movingwin, params, up_markers1, up_markers2);
        end
        m = n_used_segs{d}(n)*params.tapers(2);
        Cw8_wup{d}(n,:) = tanh(atanh(Cw8) - 1/(2*m-2));
%         Cw5_wup{d}(n,:) = tanh(atanh(Cw5) - 1/(2*m-2));
        Sw_wup{d}(n,:) = log(Sw8(:,1)) - psi(m) + log(m);
        S8_wup{d}(n,:) = log(Sw8(:,2)) - psi(m) + log(m);
%         S5_wup{d}(n,:) = log(Sw5(:,2)) - psi(m) + log(m);
        if sess_data(d).thom_elec
            Cw4_wup{d}(n,:) = tanh(atanh(Cw4) - 1/(2*m-2));
            C84_wup{d}(n,:) = tanh(atanh(C84) - 1/(2*m-2));    
            S4_wup{d}(n,:) = log(S84(:,2)) - psi(m) + log(m);
       end
       [w8_up_xcorr{d}(n,:),std_xcorr,lags] = sliding_win_xcorr(wcv_f,lf8_f,maxlag,up_markers1,up_markers2);
       if sess_data(d).thom_elec
           [w4_up_xcorr{d}(n,:),std_xcorr,lags] = sliding_win_xcorr(wcv_f,lf4_f,maxlag,up_markers1,up_markers2);
           [x84_up_xcorr{d}(n,:),std_xcorr,lags] = sliding_win_xcorr(lf8_f,lf4_f,maxlag,up_markers1,up_markers2);
       end
    end
            
    [Sw,f,varSw]=mtspectrumsegc(wcv_f,movingwin(1),params);
    [S8,f,varS8]=mtspectrumsegc(lf8_f,movingwin(1),params);
%     [S5,f,varS5]=mtspectrumsegc(lf5_f,movingwin(1),params);
    if sess_data(d).thom_elec
       [S4,f,varS4] = mtspectrumsegc(lf4_f,movingwin(1),params); 
    end
    specw(d,:) = Sw;
    spec8(d,:) = S8;
    varspecw(d,:) = varSw;
    varspec8(d,:) = varS8;
    
    nSw_up{d} = Sw_wup{d} - repmat(log(Sw)',size(Sw_wup{d},1),1);
    nSw_up{d} = nSw_up{d}./repmat(sqrt(varSw'),size(Sw_wup{d},1),1);
    nS8_up{d} = S8_wup{d} - repmat(log(S8)',size(S8_wup{d},1),1);
    nS8_up{d} = nS8_up{d}./repmat(sqrt(varS8'),size(S8_wup{d},1),1);
%     nS5_wup{d} = S5_wup{d} - repmat(log(S5)',size(S5_wup{d},1),1);
%     nS5_wup{d} = nS5_wup{d}./repmat(sqrt(varS5'),size(S5_wup{d},1),1);
    if sess_data(d).thom_elec
        nS4_up{d} = S4_wup{d} - repmat(log(S4)',size(S4_wup{d},1),1);
        nS4_up{d} = nS4_up{d}./repmat(sqrt(varS4'),size(S4_wup{d},1),1);
    end
    
    max_dur = 2;
    t = movingwin(1)/2-back_wins*movingwin(2):movingwin(2):max_dur-movingwin(1)/2;
    st = size(Cw8_wup{d},1);
    if  length(t) > st
        t(st+1:end) = [];
    end
    subplot(2,1,1)
    pcolor(t,f,Cw8_wup{d}(1:length(t),:)');shading flat
    subplot(2,1,2)
    pcolor(t,lags/Fsd,w8_up_xcorr{d}(1:length(t),:)');shading flat
    ylim([-0.2 0.2]), xl = xlim; line(xl,[0 0],'color','w')
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['G:\WC_Germany\parietal_cortical_2010\trig_mt_tfd\individ\Cw8_' cname];
    print(t_names,'-dpng')
    close

    max_dur = 2;
    t = movingwin(1)/2-back_wins*movingwin(2):movingwin(2):max_dur-movingwin(1)/2;
    st = size(Cw8_wup{d},1);
    if  length(t) > st
        t(st+1:end) = [];
    end
    subplot(2,1,1)
    pcolor(t,f,nSw_up{d}(1:length(t),:)');shading flat
    subplot(2,1,2)
    pcolor(t,f,nS8_up{d}(1:length(t),:)');shading flat
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['G:\WC_Germany\parietal_cortical_2010\trig_mt_tfd\individ\nSw8_' cname];
    print(t_names,'-dpng')
    close

    if sess_data(d).thom_elec
        max_dur = 2;
        t = movingwin(1)/2-back_wins*movingwin(2):movingwin(2):max_dur-movingwin(1)/2;
        subplot(2,1,1)
        pcolor(t,f,Cw4_wup{d}(1:length(t),:)');shading flat
        subplot(2,1,2)
        pcolor(t,lags/Fsd,w4_up_xcorr{d}(1:length(t),:)');shading flat
        ylim([-0.2 0.2]), xl = xlim; line(xl,[0 0],'color','w')
        cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
        t_names = ['G:\WC_Germany\parietal_cortical_2010\trig_mt_tfd\individ\Cw4_' cname];
        print(t_names,'-dpng')
        close
        
        max_dur = 2;
        t = movingwin(1)/2-back_wins*movingwin(2):movingwin(2):max_dur-movingwin(1)/2;
        subplot(2,1,1)
        pcolor(t,f,C84_wup{d}(1:length(t),:)');shading flat
        subplot(2,1,2)
        pcolor(t,lags/Fsd,x84_up_xcorr{d}(1:length(t),:)');shading flat
        ylim([-0.2 0.2]), xl = xlim; line(xl,[0 0],'color','w')
        cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
        t_names = ['G:\WC_Germany\parietal_cortical_2010\trig_mt_tfd\individ\84_' cname];
        print(t_names,'-dpng')
        close

            max_dur = 2;
    t = movingwin(1)/2-back_wins*movingwin(2):movingwin(2):max_dur-movingwin(1)/2;
    st = size(Cw8_wup{d},1);
    if  length(t) > st
        t(st+1:end) = [];
    end
    subplot(2,1,1)
    pcolor(t,f,nSw_up{d}(1:length(t),:)');shading flat
    subplot(2,1,2)
    pcolor(t,f,nS4_up{d}(1:length(t),:)');shading flat
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['G:\WC_Germany\parietal_cortical_2010\trig_mt_tfd\individ\nSw4_' cname];
    print(t_names,'-dpng')
    close
   
        max_dur = 2;
    t = movingwin(1)/2-back_wins*movingwin(2):movingwin(2):max_dur-movingwin(1)/2;
    st = size(Cw8_wup{d},1);
    if  length(t) > st
        t(st+1:end) = [];
    end
    subplot(2,1,1)
    pcolor(t,f,nS8_up{d}(1:length(t),:)');shading flat
    subplot(2,1,2)
    pcolor(t,f,nS4_up{d}(1:length(t),:)');shading flat
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['G:\WC_Germany\parietal_cortical_2010\trig_mt_tfd\individ\nS84_' cname];
    print(t_names,'-dpng')
    close
    end
    
end

%%
cd G:\WC_Germany\parietal_cortical_2010\
save parietal_lf8_utrig_mtgrams_200ms_xc

n = length(sess_data);
%%
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);
sup_par = intersect(parietal,superficial);
sup_deep = intersect(parietal,deep);
all = 1:n-1;

cur_set = thom_el;

%%
t_up = movingwin(1)/2-back_wins*movingwin(2):movingwin(2):max_dur-movingwin(1)/2;
max_segs = length(t_up);
used_cells = zeros(max_segs,1);
used_4cells = zeros(max_segs,1);
avg_cw8 = zeros(max_segs,length(f));
% avg_cw5 = zeros(max_segs,length(f));
avg_cw4 = zeros(max_segs,length(f));
avg_c84 = zeros(max_segs,length(f));
avg_nSw = zeros(max_segs,length(f));
avg_nS8 = zeros(max_segs,length(f));
% avg_nS5 = zeros(max_segs,length(f));
avg_nS4 = zeros(max_segs,length(f));
avg_Sw = zeros(max_segs,length(f));
avg_S8 = zeros(max_segs,length(f));
% avg_S5 = zeros(max_segs,length(f));
avg_S4 = zeros(max_segs,length(f));
% avg_xcorr = zeros(max_segs,length(lags));
for d = 1:length(cur_set)
    cur_nsegs = min(max_segs,length(n_used_segs{cur_set(d)}));
    used_cells(1:cur_nsegs) = used_cells(1:cur_nsegs)+1;
    avg_cw8(1:cur_nsegs,:) = avg_cw8(1:cur_nsegs,:) + atanh(Cw8_wup{cur_set(d)}(1:cur_nsegs,:));
%     avg_cw5(1:cur_nsegs,:) = avg_cw5(1:cur_nsegs,:) + atanh(Cw5_wup{cur_set(d)}(1:cur_nsegs,:));
    avg_nSw(1:cur_nsegs,:) = avg_nSw(1:cur_nsegs,:) + nSw_wup{cur_set(d)}(1:cur_nsegs,:);
    avg_nS8(1:cur_nsegs,:) = avg_nS8(1:cur_nsegs,:) + nS8_wup{cur_set(d)}(1:cur_nsegs,:);
%     avg_nS5(1:cur_nsegs,:) = avg_nS5(1:cur_nsegs,:) + nS5_wup{cur_set(d)}(1:cur_nsegs,:);
    avg_Sw(1:cur_nsegs,:) = avg_Sw(1:cur_nsegs,:) + Sw_wup{cur_set(d)}(1:cur_nsegs,:);
    avg_S8(1:cur_nsegs,:) = avg_S8(1:cur_nsegs,:) + S8_wup{cur_set(d)}(1:cur_nsegs,:);
%     avg_S5(1:cur_nsegs,:) = avg_S5(1:cur_nsegs,:) + S5_wup{cur_set(d)}(1:cur_nsegs,:);
    if sess_data(cur_set(d)).thom_elec
        used_4cells(1:cur_nsegs) = used_4cells(1:cur_nsegs)+1;
        avg_nS4(1:cur_nsegs,:) = avg_nS4(1:cur_nsegs,:) + nS4_wup{cur_set(d)}(1:cur_nsegs,:);
        avg_S4(1:cur_nsegs,:) = avg_S4(1:cur_nsegs,:) + S4_wup{cur_set(d)}(1:cur_nsegs,:);
        avg_cw4(1:cur_nsegs,:) = avg_cw4(1:cur_nsegs,:) + atanh(Cw4_wup{cur_set(d)}(1:cur_nsegs,:));
        avg_c84(1:cur_nsegs,:) = avg_c84(1:cur_nsegs,:) + atanh(C84_wup{cur_set(d)}(1:cur_nsegs,:));
    end
%     avg_xcorr(1:cur_nsegs,:) = avg_xcorr(1:cur_nsegs,:) + w8_wup_xcorr{cur_set(d)}(1:cur_nsegs,:);
end

avg_cw8 = tanh(avg_cw8./repmat(used_cells,1,length(f)));
avg_cw8(used_cells < 5,:) = nan;
% avg_cw5 = tanh(avg_cw5./repmat(used_cells,1,length(f)));
% avg_cw5(used_cells < 5,:) = nan;
avg_nSw = avg_nSw./repmat(used_cells,1,length(f));
avg_nSw(used_cells < 5,:) = nan;
avg_nS8 = avg_nS8./repmat(used_cells,1,length(f));
avg_nS8(used_cells < 5,:) = nan;
% avg_nS5 = avg_nS5./repmat(used_cells,1,length(f));
% avg_nS5(used_cells < 5,:) = nan;
avg_Sw = avg_Sw./repmat(used_cells,1,length(f));
avg_Sw(used_cells < 5,:) = nan;
avg_S8 = avg_S8./repmat(used_cells,1,length(f));
avg_S8(used_cells < 5,:) = nan;
% avg_S5 = avg_S5./repmat(used_cells,1,length(f));
% avg_S5(used_cells < 5,:) = nan;
avg_cw4 = tanh(avg_cw4./repmat(used_4cells,1,length(f)));
avg_cw4(used_4cells < 5,:) = nan;
avg_c84 = tanh(avg_c84./repmat(used_4cells,1,length(f)));
avg_c84(used_4cells < 5,:) = nan;
avg_nS4 = avg_nS4./repmat(used_4cells,1,length(f));
avg_nS4(used_4cells < 5,:) = nan;
avg_S4 = avg_S4./repmat(used_4cells,1,length(f));
avg_S4(used_4cells < 5,:) = nan;

%%
figure
pcolor(t_up,f,avg_c84');shading flat
caxis([-0.1 0.7])
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Up-trig','fontsize',18)
colorbar
%%
figure
pcolor(t_up,f,avg_cw8');shading flat
caxis([-0.1 0.2])
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Up-trig','fontsize',18)
colorbar
%%
figure
pcolor(t_up,f,avg_cw4');shading flat
caxis([-0.1 0.2])
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Up-trig','fontsize',18)
colorbar

%%
figure
pcolor(t_up,f,avg_nS4');shading flat
caxis([-1. 0.7])
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Up-trig','fontsize',18)
colorbar
%%
figure
pcolor(t_up,f,avg_S4');shading flat
% caxis([-1. 0.7])
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Up-trig','fontsize',18)
colorbar
%%
figure
pcolor(t_up,f,avg_S8');shading flat
caxis([-15 -11])
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Up-trig','fontsize',18)
colorbar

%%
figure
pcolor(t_up,f,avg_Sw');shading flat
caxis([-9 -4.5])
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Up-trig','fontsize',18)
colorbar

%%
% figure
% subplot(2,1,1)
pcolor(t_up,f,avg_nSw');shading flat
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('MP','fontsize',18)
colorbar

% subplot(2,1,2)
figure
pcolor(t_up,f,avg_nS8');shading flat
xlabel('Time (s)','fontsize',16)
ylabel('Lag (s)','fontsize',16)
title('LFP','fontsize',18)
colorbar
