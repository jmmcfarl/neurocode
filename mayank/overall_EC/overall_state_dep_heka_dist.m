clear all

cd G:\WC_Germany\overall_EC\
addpath('G:\WC_Germany\hsmm_state_detection\')
load G:\WC_Germany\overall_EC\overall_EC_dir.mat
load heka_alignment.mat
dsf = 8;
Fsd = 2016/dsf;
% fs_desired = 500;

amp_range = linspace(-85,-20,500);

%loop over just the MEC cells until get the full LEC heka data
for d = 64:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    if ~isempty(sess_data(d).heka_dir)
        load used_data wcv_minus_spike
        load ec_hmm_state_seq
        load aligned_heka

        mp_state_seq = hmm_bbstate_seq;

        wcv_d = downsample(wcv_minus_spike,dsf);
        dc_data = downsample(dc_data,dsf);
        dc_data = dc_data(:);
        if mean(abs(dc_data)) < 1
            dc_data = dc_data*100;
        end
        dc_time = downsample(dc_time,dsf);
        ac_time = downsample(ac_time,dsf);

        [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(wcv_d));

        up_dc_sig = [];
        down_dc_sig = [];
        for i = 1:hmm.Nsegs
            cur_up_trans{i} = find(mp_state_seq{i}(1:end-1) == 1 & mp_state_seq{i}(2:end) == 2);
            cur_down_trans{i} = find(mp_state_seq{i}(1:end-1) == 2 & mp_state_seq{i}(2:end) == 1);
            cur_up_trans{i}(cur_up_trans{i} > length(dc_time)) = [];
            cur_down_trans{i}(cur_down_trans{i} > length(dc_time)) = [];
            cur_up_trans{i}(cur_up_trans{i} > length(ac_time)) = [];
            cur_down_trans{i}(cur_down_trans{i} > length(ac_time)) = [];            
            cur_up_trans{i}(cur_up_trans{i} > cur_down_trans{i}(end)) = [];
            cur_down_trans{i}(cur_down_trans{i} < cur_up_trans{i}(1)) = [];
            cur_up_trans{i} = cur_up_trans{i} + new_seg_inds(i,1)-1;
            cur_down_trans{i} = cur_down_trans{i} + new_seg_inds(i,1)-1;
        
            dc_up_ind{i} = nan(size(cur_up_trans{i}));
            dc_down_ind{i} = nan(size(cur_up_trans{i}));
            for s = 1:length(cur_up_trans{i})
                cur_up_ind = find(dc_time > ac_time(cur_up_trans{i}(s)),1,'first');
                cur_down_ind = find(dc_time > ac_time(cur_down_trans{i}(s)),1,'first');
                if ~isempty(cur_up_ind) && ~isempty(cur_down_ind)
                    dc_up_ind{i}(s) = cur_up_ind;
                    dc_down_ind{i}(s) = cur_down_ind;
                    up_dc_sig = [up_dc_sig; dc_data(dc_up_ind{i}(s):dc_down_ind{i}(s))];
                end
            end
            for s = 1:length(cur_up_trans{i})-1
                if ~isnan(dc_up_ind{i}(s+1)) && ~isnan(dc_down_ind{i}(s))
                    down_dc_sig = [down_dc_sig; dc_data(dc_down_ind{i}(s):dc_up_ind{i}(s+1))];
                end
            end

        end

        upstate_heka_dens(d,:) = ksdensity(up_dc_sig,amp_range);
        downstate_heka_dens(d,:) = ksdensity(down_dc_sig,amp_range);
        upstate_mean(d) = mean(up_dc_sig);
        downstate_mean(d) = mean(down_dc_sig);
        upstate_var(d) = var(up_dc_sig);
        downstate_var(d) = var(down_dc_sig);

        plot(amp_range,upstate_heka_dens(d,:),'r'), hold on
        plot(amp_range,downstate_heka_dens(d,:))
        t_names = ['G:\WC_Germany\overall_EC\state_dep_heka_dist\' s_name];
        print(t_names,'-dpng')
        close
    end
    clear ac_time* dc_*
    
end

cd G:\WC_Germany\overall_EC
save heka_state_dep_data upstate* downstate*