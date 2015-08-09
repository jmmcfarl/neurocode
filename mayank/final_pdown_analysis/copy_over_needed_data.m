clear all
close all

addpath('~/James_scripts/mayank/final_pdown_analysis/');

database_file = '/Users/james/Analysis/Mayank/final_pdown_analysis/compiled_data.mat';
% database_file = '/Users/james/Analysis/Mayank/final_pdown_analysis/compiled_corticalMP_data.mat';
load(database_file);
target_base_dir = '/Users/james/Analysis/Mayank/';

% target_files = {'used_data.mat','pa_hsmm_state_seq7_combined_fin_nd.mat',...
%     'pa_hsmm_state_seq_combined_fin_nd.mat','mua_data3.mat',...
%     'sync_times.mat','spike_time_jmm.mat'};

target_files = {};
% target_files = {'pa_hsmm_state_seq7_combined_fin_nd.mat',...
%     'pa_hsmm_state_seq_combined_fin_nd.mat','allEC_ctx_period_data_hsmm.mat'};
% target_files = {'allEC_ctx_period_data_hsmm.mat'};
    %%
for dd = 1:length(data)
    
    cur_dir = data(dd).dir;
    
    new_dir = map_to_new_drive_locs(cur_dir);
    
    cd(new_dir)
    pwd
    
    out_dir = strrep(cur_dir,'C:\',target_base_dir);
    out_dir = strrep(out_dir,'\','/');
    if ~exist(out_dir,'dir');
        system(sprintf('mkdir -p %s',out_dir));
    end
    
    for ii = 1:length(target_files)
        dest_file = sprintf('%s/%s',out_dir,target_files{ii});
        if exist(target_files{ii},'file') && ~exist(dest_file,'file')
            from_name = sprintf('./%s',target_files{ii});
            fprintf('Copying %s to %s\n',from_name,out_dir);
            system(sprintf('cp %s %s',from_name,out_dir));
        elseif ~exist(target_files{ii},'file')
            warning('%s not found',target_files{ii});
        end
    end
    
    data(dd).new_dir = out_dir;
end

%%
save(database_file,'data','no_cell','unclear_uds','clear*')