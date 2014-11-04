clear all
close all
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd_dist.mat
addpath('C:\WC_Germany\persistent_2010\')
%%
% true_heka_fs = 1.999834431847603e+004;
true_heka_fs = 1.99984e+004;
% true_heka_fs = 2e+004;

for d = 77:81
    disp(combined_heka{d});
    cd(combined_dir{d})
    load ./used_data wcv
    load ./sync_times
    hloc = combined_heka{d};
    temp_data = load(hloc);
    hloc(1:find(hloc=='\',1,'last')) = [];
    if isfield(temp_data,'data')
        dc_data = temp_data.data;
    elseif isfield(temp_data,'dataVec')
        dc_data = temp_data.dataVec;
    else
        eval(['dc_data = temp_data.' hloc '_MP;']);
    end
    if isfield(temp_data,'time')
        dc_times = temp_data.time;
    elseif isfield(temp_data,'timeVec')
        dc_times = temp_data.timeVec;
    elseif isfield(temp_data,strcat(hloc,'_sampletimes'))
        eval(['dc_times = temp_data.' hloc '_sampletimes;']);
    else
        dc_times = (1:length(dc_data))/true_heka_fs;
    end
    clear temp_data
    if strcmp(combined_heka_type{d},'sweep')
        [sweep_offsets{d},sweep_maxcorrs{d},sweep_starts{d},sweep_stops{d}] = align_dc_ac_sigs_v2(dc_data,dc_times,wcv);
        ac_time = (1:length(wcv))/2016;
        dc_time = dc_times;
        for i = 1:length(sweep_starts{d})
            dc_time(sweep_starts{d}(i):sweep_stops{d}(i)) = dc_time(sweep_starts{d}(i):sweep_stops{d}(i))+sweep_offsets{d}(i)-dc_time(sweep_starts{d}(i));
        end
        save aligned_heka dc_data dc_time ac_time
    elseif strcmp(combined_heka_type{d},'cont')
        %         interp_ax = synct(1):median(diff(synct)):synct(end);
        %         wcv_int = interp1(synct,wcv,interp_ax);
        wcv_int = wcv;
        [dc_offset(d),dc_maxcorr(d)] = align_dc_ac_sigs_initial_v2(dc_data,dc_times,wcv_int);
        dc_time = dc_times + dc_offset(d);
        ac_time = (1:length(wcv_int))/2016;
        %         ac_time = interp1(interp_ax,ac_time,synct);
        if min(diff(dc_time)) <= 0
            error('Alignment Problem!')
        end
        save aligned_heka dc_data dc_time ac_time
    else
        error('invalid heka type')
    end
end

% cd C:\WC_Germany\sven_thomas_combined\
% save combined_heka_alignment sweep_* dc_offset dc_maxcorr
%%
% for d = 1:51
%     d
%         cd(combined_dir{d})
%     load ./used_data wcv
%     load ./aligned_heka
%
%     t = (1:length(wcv))/2016;
%     plot(t,zscore(wcv))
%     hold on
%
%     plot(downsample(dc_time,10),zscore(downsample(dc_data,10)),'r')
%     pause
%     close all
% end