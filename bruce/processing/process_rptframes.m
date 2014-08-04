clear all
Expt_name = 'M296';
Expt_num = str2num(Expt_name(2:end));
if Expt_num > 280 
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end
cd(data_dir);

if Expt_name(1) == 'M'
    rec_type = 'LP';
elseif Expt_name(1) == 'G'
    rec_type = 'UA';
end

if strcmp(rec_type,'LP')
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
elseif strcmp(rec_type,'UA')
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
end

%%
% tot_blocks = length(Expts);
% cd([data_dir '/stims/']);
% for bb = 1:tot_blocks
%    cur_fname = sprintf('Expt%d_stim.mat',bb);
%    if exist(cur_fname,'file')
%        load(cur_fname);
%        new_name = sprintf('Expt%d_oldStim.mat',bb);
%        save(new_name,'left_stim_mats','right_stim_mats','stim_nframes','stim_binoc');
%        
%        fprintf('Loaded %s, saved as %s\n',cur_fname,new_name);
%        clear left_stim_mats right_stim_mats stim_nframes stim_binoc
%    end
%     
% end
% 
%%

if strcmp(rec_type,'LP')
    switch Expt_num
        case 266
            bar_ori = 80;
        case 270
            bar_ori = 60;
        case 275
            bar_ori = 135;
        case 277
            bar_ori = 70;
        case 281
            bar_ori = 140;
        case 287
            bar_ori = 90;
        case 289
            bar_ori = 160;
        case 294
            bar_ori = 40;
        case 296
            bar_ori = 45;
    end
end


load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

%%
stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;


ignore_blocks = [];

%% SELECT BLOCKS FOR ANALYSIS
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi','rls.AllSacB'};
expt_names = cell(1,length(Expts));
expt_dds = nan(1,length(Expts));
expt_bar_ori = nan(1,length(Expts));
expt_sac_dir = nan(1,length(Expts));
expt_Fr = nan(1,length(Expts));
expt_imback = nan(1,length(Expts));
included_type = false(1,length(Expts));
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;

if strcmp(rec_type,'LP')
    expt_bar_ori(expt_bar_ori > 360) = bar_ori;
end

cur_block_set = find(included_type & expt_Fr == 1 & expt_bar_ori == bar_ori);

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];
n_blocks = length(cur_block_set);


%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_stim_mat = [];
all_t_axis = [];
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:n_blocks;
    fprintf('Expt %d Block %d of %d\n',Expt_num,ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= 0);
    end
            
    fname = sprintf('%s/stims/Expt%d_oldStim',data_dir,cur_block);
%     fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        if isfield(Expts{cur_block}.Trials(tt),'rptframes')
            cur_rpt_frames = Expts{cur_block}.Trials(tt).rptframes;
            
%             cur_rpt_frames(cur_rpt_frames == 1) = [];
            
            if ~isempty(cur_rpt_frames)
                n_frames = size(left_stim_mats{use_trials(tt)},1);
                frame_vec = 1:n_frames;
                
                for ii = length(cur_rpt_frames):-1:1
                    frame_vec = [frame_vec(1:cur_rpt_frames(ii)) cur_rpt_frames(ii) frame_vec(cur_rpt_frames(ii)+1:end)];
                end
                
                left_stim_mats{tt} = left_stim_mats{tt}(frame_vec,:);
                right_stim_mats{tt} = right_stim_mats{tt}(frame_vec,:);
                              
                fprintf('Blk %d, Trial %d, replaced %d rptframes\n',cur_block,tt,length(cur_rpt_frames));
            end
        end
    end
    
    out_fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    save(out_fname,'left_stim_mats','right_stim_mats','stim_nframes','stim_binoc');
    
end









