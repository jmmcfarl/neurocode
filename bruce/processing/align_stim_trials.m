clear all

expt_name = 'M297';
%%

% cd(['~/Data/bruce/' expt_name]);
% cd(['/media/NTlab_data2/Data/bruce/' expt_name]);
cd(['/media/NTlab_data3/Data/bruce/' expt_name]);
if expt_name(1) == 'G'
    load(sprintf('jbe%sExpts.mat',expt_name));
elseif expt_name(1) == 'M'
    load(sprintf('lem%sExpts.mat',expt_name));
end
cd stims/
load ./stim_data.mat

%%
exclude_expt_types = {'rls.orXme','rls.orRC'};
%%

% trial_is_used = false(size(trial_trial_ids));

n_expts = length(Expts);
expt_binoc = nan(n_expts,1);
expt_has_ds = nan(n_expts,1);
expt_npix = nan(n_expts,1);
expt_ntrials = nan(n_expts,1);
expt_trial_nmismatch = nan(n_expts,1);
expt_stimmatches_found = cell(n_expts,1);
expt_rc_matches = cell(n_expts,1);
for ii = 1:n_expts
    % for ii = 17
    if ~isempty(Expts{ii})
        if ~strcmp(Expts{ii}.Header.expname,exclude_expt_types)
            
            fprintf('Aligning expt %d of %d\n',ii,n_expts);
            cur_trial_ids = [Expts{ii}.Trials(:).id];
            cur_trial_durs = [Expts{ii}.Trials(:).dur];
            expt_first_trial(ii) = cur_trial_ids(1);
            expt_last_trial(ii) = cur_trial_ids(end);
            
            n_trials = length(cur_trial_ids);
            n_un_trials = length(unique(cur_trial_ids));
            fprintf('%d trials, %d unique\n',n_trials,n_un_trials);
            cur_trial_ids = unique(cur_trial_ids);
            if length(cur_trial_ids) < 5
                fprintf('Too few trials, skipping block\n');
                cur_trial_ids = [];
                n_un_trials = 0;
            end
            
            stim_index = nan(n_un_trials,1);
            left_stim_mats = cell(n_un_trials,1);
            right_stim_mats = cell(n_un_trials,1);
            expt_seed_nums = nan(n_un_trials,1);
            stim_nframes = nan(n_un_trials,1);
            stim_binoc = nan(n_un_trials,1);
            expt_rc_matches{ii} = nan(n_un_trials,1);
            %         expt_trial_ds{ii} = nan(n_un_trials,1);
            for jj = 1:n_un_trials
                %find rc trials with matching trial ids
                match = find(trial_trial_ids == cur_trial_ids(jj));
                
                if strcmp(expt_name,'G086')
                    match(trial_file_ids(match) == 23) = []; %rls.rc193 just has two trials and creates problems
                end
                
                %             match(trial_is_used(match)) = [];
                cur_n_matches(jj) = length(match);
                if length(match) > 1
                    %                 fprintf('Warning multiple trial id matches\n');
                    %take first match that hasn't already been used
                    if ii < 14
                        match = match(1);
                    else
                        match = match(2);
                    end
                end
                if ~isempty(match)
                    stim_index(jj) = match;
                    %                 trial_is_used(match) = true; %mark trial as used in stim files
                    cur_file_id = trial_file_ids(match);
                    cur_frame_set = find(frame_file_ids == cur_file_id);
                    cur_frame_inds = find(frame_trial_nums(cur_frame_set) == cur_trial_ids(jj));
                    
                    cur_left_pix = all_left_pix{cur_file_id}(cur_frame_inds,:);
                    cur_right_pix = all_right_pix{cur_file_id}(cur_frame_inds,:);
                    
                    expt_seed_nums(jj) = trial_seed_nums(match);
                    expt_rc_matches{ii}(jj) = cur_file_id;
                    %                 expt_trial_ds{ii}(jj) = trial_is_ds(match);
                    left_stim_mats{jj} = cur_left_pix;
                    right_stim_mats{jj} = cur_right_pix;
                    if all(cur_left_pix(:) == cur_right_pix(:))
                        stim_binoc(jj) = 0;
                    else
                        stim_binoc(jj) = 1;
                    end
                    stim_nframes(jj) = length(cur_frame_inds);
                else
                    %             fprintf('NO MATCH FOUND FOR THIS TRIAL!\n');
                end
            end
            expt_rc_index{ii} = stim_index;
            expt_trial_nmismatch(ii) = sum(isnan(expt_rc_matches{ii}));
            expt_ntrials(ii) = n_un_trials;
            
            fprintf('%d of %d trials not matched\n',expt_trial_nmismatch(ii),expt_ntrials(ii));
            if expt_ntrials(ii)-expt_trial_nmismatch(ii) > 0
                
                backwards_steps = find(expt_rc_matches{ii}(1:end-1) < expt_rc_matches{ii}(2:end));
                n_backsteps = length(backwards_steps);
                if n_backsteps > 0
                    fprintf('%d backwards matches found\n',n_backsteps);
                    expt_rc_matches{ii}(backwards_steps+1) = nan;
                    stim_nframes(backwards_steps+1) = nan;
                    stim_binoc(backwards_steps+1) = nan;
                    for jjj = 1:n_backsteps
                        left_stim_mats{backwards_steps(jjj)+1} = nan;
                        right_stim_mats{backwards_steps(jjj)+1} = nan;
                    end
                end
                
                rc_files_linked = unique(expt_rc_matches{ii}(~isnan(expt_rc_matches{ii})));
                %         linked_rc_names = cell2mat(all_file_names(rc_files_linked));
                linked_rc_names = all_file_names(rc_files_linked);
                fprintf('Linked expt %d with rc files:  ',ii);
                for cc = 1:length(linked_rc_names)
                    fprintf('%s  ',linked_rc_names{cc});
                end
                fprintf('\n');
                expt_npix(ii) = size(cur_left_pix,2);
                
                expt_binoc(ii) = max(stim_binoc);
                if expt_binoc(ii) ~= any(rc_is_binoc(rc_files_linked))
                    error('Problem with linking, binoc mismatch');
                end
                
                un_ds = unique(rc_has_ds(rc_files_linked));
                if length(un_ds) > 1
                    expt_has_ds(ii) = -1;
                else
                    expt_has_ds(ii) = un_ds;
                end
                if expt_has_ds(ii)
                    fprintf('Expt %d has Guided Sac\n',ii);
                else
                    fprintf('Expt %d does NOT have Guided Sac\n',ii);
                end
                
            end
            
            %for M296, the rpt trials got stims set to nan, so manually set
            %them to their true values (stored in rpt_stim.mat)
            if str2num(expt_name(2:end)) == 296
               load ./rpt_stim.mat
               rpt_seed = find(expt_seed_nums == 1803);
                
               for jj = 1:length(rpt_seed)
                  cur_tlen = size(left_stim_mats{rpt_seed(jj)},1); 
                  left_stim_mats{rpt_seed(jj)} = -rpt_stim(1:cur_tlen,:);
                  right_stim_mats{rpt_seed(jj)} = -rpt_stim(1:cur_tlen,:); %contrast was -1 so the positive and negative values were flipped
               end
               
            end
            
            %%
            fprintf('\n');
            sname = sprintf('Expt%d_stim',ii);
            save(sname,'stim_nframes','stim_binoc','left_stim_mats','right_stim_mats');
            
        end
    end
end

%% SAVE SOME META-DATA ABOUT THE EXPTS STIMULI AND MATCHING
save expt_data expt*