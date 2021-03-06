clear all
close all

addpath('~/James_scripts/mayank/overall_EC/');
addpath('~/James_scripts/mayank/final_pdown_analysis/');

% load C:/WC_Germany/final_pdown_analysis/compiled_data.mat
load ~/Analysis/Mayank/final_pdown_analysis/compiled_data.mat

true_heka_fs = 1.99984e+004;
% spk_win_t = .004;
for d = 81:length(data)
    if ~isempty(data(d).heka_data)
        cur_dir = data(d).dir;
        new_dir = map_to_new_drive_locs(cur_dir);
        cd(new_dir)
        pwd
        
        load ./used_data wcv
        load ./sync_times.mat
        
        %         heka_data = load(data(d).heka_data);
        new_heka_loc = map_to_new_drive_locs(data(d).heka_data,1);
        if ~exist(new_heka_loc,'file')
            new_heka_loc = map_to_new_drive_locs(data(d).heka_data,1,1);
        end
         heka_data = load(new_heka_loc);
        dat_fields = fieldnames(heka_data);
        clear MP_loc time_loc
        for ii = 1:length(dat_fields)
            if strcmp(dat_fields{ii}(end-1:end),'MP')
                MP_loc = ii;
            elseif length(dat_fields{ii}) > 4 && strcmp(dat_fields{ii}(end-4:end),'times')
                time_loc = ii;
            end
        end
        if isfield(heka_data,'data')
            dc_data = heka_data.data;
        elseif isfield(heka_data,'dataVec')
            dc_data = heka_data.dataVec;
        elseif isfield(heka_data,'dc_data')
            dc_data = heka_data.dc_data;
        else
            dc_data = heka_data.(dat_fields{MP_loc});
        end
        if isfield(heka_data,'time')
            dc_times = heka_data.time;
        elseif isfield(heka_data,'timeVec')
            dc_times = heka_data.timeVec;
        elseif length(dat_fields) > 1
            temp = fieldnames(heka_data);
            dc_times = heka_data.(temp{time_loc});
        else
            dc_times = (1:length(dc_data))/true_heka_fs;
        end
        clear heka_data
        
        %     load ./spike_time_jmm.mat
        
        if strcmp(data(d).heka_type,'sweep')
            [sweep_offsets,sweep_maxcorrs,sweep_starts,sweep_stops] = align_dc_ac_sigs_v2(dc_data,dc_times,wcv);
            ac_time = (1:length(wcv))/2016;
            dc_time = dc_times;
            dc_dt = median(diff(dc_time));
            for i = 1:length(sweep_starts)
                dc_time(sweep_starts(i):sweep_stops(i)) = dc_time(sweep_starts(i):sweep_stops(i))+sweep_offsets(i)-dc_time(sweep_starts(i));
            end
            
            %         spike_times = ac_time(spkid)-0.001;
            %
            %         dc_spike_inds = round(interp1(dc_time,1:length(dc_time),spike_times));
            %         spike_times(isnan(dc_spike_inds)) = [];
            %         dc_spike_inds(isnan(dc_spike_inds)) = [];
            %
            %         spike_error = abs(dc_time(dc_spike_inds) - spike_times');
            %         bad = find(spike_error > 2/dc_dt);
            %         dc_spike_inds(bad) = [];
            %         spk_win = round(spk_win_t/dc_dt);
            %         interp_inds = bsxfun(@plus,dc_spike_inds',1:spk_win);
            %         interp_inds = interp_inds(:);
            %         used_inds = setdiff(1:length(dc_data),interp_inds);
            %         dc_data_spksub = interp1(used_inds,dc_data(used_inds),1:length(dc_data));
            
            out_file_name = [data(d).new_dir '/aligned_heka'];
            save(out_file_name,'dc_data','dc_time','ac_time','sweep*');
%             save aligned_heka dc_data dc_time ac_time
            fprintf('Aligned with avg correlation %.3f\n',mean(sweep_maxcorrs));
        elseif strcmp(data(d).heka_type,'cont')
            wcv_int = wcv;
            [dc_offset,dc_maxcorr] = align_dc_ac_sigs_initial_v2(dc_data,dc_times,wcv_int);
            dc_time = dc_times + dc_offset;
            ac_time = (1:length(wcv_int))/2016;
            if min(diff(dc_time)) <= 0
                error('Alignment Problem!')
            end
            
            %         spike_times = ac_time(spkid)-0.001;
            %
            %         dc_spike_inds = round(interp1(dc_time,1:length(dc_time),spike_times));
            %         spike_times(isnan(dc_spike_inds)) = [];
            %         dc_spike_inds(isnan(dc_spike_inds)) = [];
            %
            %         spike_error = abs(dc_time(dc_spike_inds) - spike_times);
            %         bad = find(spike_error > 2/dc_dt);
            %         dc_spike_inds(bad) = [];
            %         spk_win = round(spk_win_t/dc_dt);
            %         interp_inds = bsxfun(@plus,dc_spike_inds',1:spk_win);
            %         interp_inds = interp_inds(:);
            %         used_inds = setdiff(1:length(dc_data),interp_inds);
            %         dc_data_spksub = interp1(used_inds,dc_data(used_inds),1:length(dc_data));
            
            out_file_name = [data(d).new_dir '/aligned_heka'];
            save(out_file_name,'dc_data','dc_time','ac_time','dc_offset','dc_maxcorr');
%             save aligned_heka dc_data dc_time ac_time
            fprintf('Aligned with correlation %.3f\n',dc_maxcorr);
        else
            error('invalid heka type')
        end
        
        clear dc_data  dc_times
%     else
%         dc_offset(d) = nan;
%         dc_maxcorr(d) = nan;
    end
end

%%
% cd C:\WC_Germany\final_pdown_analysis\
% cd ~/Analysis/Mayank/final_pdown_analysis/
% save overall_heka_alignment_stats dc_offset dc_maxcorr sweep_offsets sweep_maxcorrs sweep_starts sweep_stops
%%
