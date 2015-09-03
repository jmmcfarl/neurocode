clear all
close all
load C:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('C:\WC_Germany\persistent_2010\')
%%
% true_heka_fs = 1.999834431847603e+004;
true_heka_fs = 1.99984e+004;
% true_heka_fs = 2e+004;

used_data = [sup_mec sup_lec];

% lec = find_struct_field_vals(sess_data,'region','LEC');
for d = used_data
    disp(d);
    cdir = sess_data(d).directory;
    cdir(1) = 'C';
    cd(cdir);
    if exist('./used_data.mat','file')
        load ./used_data wcv
        hloc = sess_data(d).heka_dir;
        %     if ~exist('./aligned_heka.mat','file')
        if ~isempty(hloc)
            hloc(1) = 'C';
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
            if strcmp(sess_data(d).heka_type,'sweep')
                %                 [best_match{d},match_qual,sweep_starts{d},sweep_stops{d}] = align_dc_ac_sigs(dc_data,dc_times,wcv_minus_spike);
                [sweep_offsets{d},sweep_maxcorrs{d},sweep_starts{d},sweep_stops{d}] = align_dc_ac_sigs_v2(dc_data,dc_times,wcv);
                ac_time = (1:length(wcv))/2016;
                dc_time = dc_times;
                for i = 1:length(sweep_starts{d})
                    dc_time(sweep_starts{d}(i):sweep_stops{d}(i)) = dc_time(sweep_starts{d}(i):sweep_stops{d}(i))+sweep_offsets{d}(i)-dc_time(sweep_starts{d}(i));
                end
                %                 dc_time = nan(size(dc_data));
                %                 best_match{d}(best_match{d} > length(ac_time)) = length(ac_time);
                %                 for i = 1:size(best_match{d},2)
                %                     cdc_inds = sweep_starts{d}(i):sweep_stops{d}(i);
                %                     dc_time(cdc_inds) = linspace(ac_time(best_match{d}(1,i)),...
                %                         ac_time(best_match{d}(2,i)),length(cdc_inds));
                %                 end
                %                 if min(diff(dc_time)) <= 0
                %                     error('Alignment Problem!')
                %                 end
                save aligned_heka dc_data dc_time ac_time
            elseif strcmp(sess_data(d).heka_type,'cont')
                %                 [ac_offset(d)] = align_dc_ac_sigs_initial(dc_data,dc_times,wcv_minus_spike);
                [dc_offset(d),dc_maxcorr(d)] = align_dc_ac_sigs_initial_v2(dc_data,dc_times,wcv);
                dc_time = dc_times + dc_offset(d);
                ac_time = (1:length(wcv))/2016;
                if min(diff(dc_time)) <= 0
                    error('Alignment Problem!')
                end
                save aligned_heka dc_data dc_time ac_time
            else
                error('invalid heka type')
            end
        end
    end
end

% cd F:\WC_Germany\Overall_EC\
% save heka_alignment sweep_* dc_offset dc_maxcorr