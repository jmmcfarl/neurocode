clear all

load C:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')
addpath('C:\WC_Germany\persistent_2010\')
addpath('C:\Code\smoothing\software\')
addpath('C:\Code\general\')

% lec = find_struct_field_vals(sess_data,'region','LEC');
% lec(2) = [];
% used_data = lec;

used_data = [sup_lec sup_mec];
sess_data = sess_data(used_data);

spike_thresh = 1e4;
spike_rem_width = 3; %in ms

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'C';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    %     if ~exist('heka_aligned_spksub.mat','file')
    if ~isempty(sess_data(d).heka_dir)
        load ./aligned_heka
        if mean(abs(dc_data)) < 1
            dc_data = dc_data*100;
        end
        
        heka_spk_locs = [];
        if strcmp(sess_data(d).heka_type,'sweep')
            dc_data = dc_data(:);
            dtimes = diff(dc_time(:));
            dt = median(dtimes);
            dc_Fs = 1/dt;
            spike_rem_samps = round(spike_rem_width*1e-3*dc_Fs);
            dtimes = [Inf; dtimes];
            sweep_start_inds = find(dtimes > 0.1);
            sweep_stop_inds = [sweep_start_inds(2:end)-1; length(dc_time)];
            
            for i = 1:length(sweep_start_inds)
                cur_deriv = [0; diff(dc_data(sweep_start_inds(i):sweep_stop_inds(i)))]/dt;
                cur_spikes = find(cur_deriv(1:end-1) < spike_thresh & cur_deriv(2:end) > spike_thresh);
                cur_spk_inds = cur_spikes(:) + sweep_start_inds(i)-1;
                for j = 1:spike_rem_samps-1
                    heka_spk_locs = [heka_spk_locs; cur_spk_inds+j];
                end
                heka_spk_locs(heka_spk_locs > sweep_stop_inds(i)) = [];
            end
            dc_samps = 1:length(dc_data);
            dc_nospk_samps = dc_samps; dc_nospk_samps(heka_spk_locs) = [];
            dc_data_spksub = interp1(dc_nospk_samps,dc_data(dc_nospk_samps),dc_samps);
            
        elseif strcmp(sess_data(d).heka_type,'cont')
            dt = median(diff(dc_time));
            dc_Fs = 1/dt;
            cur_deriv = [0; diff(dc_data(:))]/dt;
            cur_spikes = find(cur_deriv(1:end-1) < spike_thresh & cur_deriv(2:end) > spike_thresh);
            for j = 1:spike_rem_samps-1
                heka_spk_locs = [heka_spk_locs; cur_spk_inds+j];
            end
            heka_spk_locs(heka_spk_locs > length(dc_data)) = [];
            
            dc_samps = 1:length(dc_data);
            dc_nospk_samps = dc_samps; dc_nospk_samps(heka_spk_locs) = [];
            dc_data_spksub = interp1(dc_nospk_samps,dc_data(dc_nospk_samps),dc_samps);
            
        else
            error('No Heka data');
        end
        
        save heka_aligned_spksub dc_data_spksub
    end
    %     end
end