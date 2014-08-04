clear all
close all

addpath('~/James_scripts/bruce/eye_tracking_finalcode/');
% expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','G085','G086','G087','G088','G089','G091','G093','G095'};
expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
print_out = false;

fig_dir = '/home/james/Analysis/bruce/ET_final/';

% use_ori = [true true];
use_ori = [true false];
% use_ori = [false true];

n_expts = length(expt_list);

%%
all_su_data = [];
for ee = 1:n_expts
    
    cur_expt_name = expt_list{ee};
    cur_expt_num = str2num(cur_expt_name(2:end));
    
    fprintf('Analyzing %s, expt %d of %d\n',cur_expt_name,ee,n_expts);
    if strcmp(cur_expt_name,'M270')
        scale_fac = 1.72;
    else
        scale_fac = 1;
    end
    
    if cur_expt_num < 280
    data_dir = ['~/Data/bruce/' cur_expt_name];
    elseif cur_expt_num < 289
    data_dir = ['/media/NTlab_data2/Data/bruce/' cur_expt_name];
    else
    data_dir = ['/media/NTlab_data3/Data/bruce/' cur_expt_name];
    end
    if cur_expt_name(1) == 'G'
        load([data_dir sprintf('/jbe%sExpts.mat',cur_expt_name)]); %load in Expts struct
    else
        load([data_dir sprintf('/lem%sExpts.mat',cur_expt_name)]); %load in Expts struct
    end
    uu = find(cellfun(@(X) length(X),Expts) > 0,1);
    rf_pos = Expts{uu}.Stimvals.rf(1:2)/scale_fac;
       
   if cur_expt_name(1) == 'G'       
       anal_dir = ['~/Analysis/bruce/' cur_expt_name '/ET_final_imp/'];
       cd(anal_dir);
       
       clear dit_*
       %for horizontal
%        mod_data_name = 'monoc_eyecorr_hbar_mods.mat';
%        anal_name = 'monoc_eyecorr_hbar2.mat';
       mod_data_name = 'full_eyetrack_initmods.mat';
       anal_name = 'full_eyetrack.mat';
       if exist(anal_name,'file') && use_ori(1)
           hor_data = load(anal_name,'it_mods','it_R2','dit_mods*','dit_R2*','et_tr_set');
           hor_mod_data = load(mod_data_name,'all_mod_SU*');
       else
           hor_mod_data.all_mod_SUnum = [];
           hor_data.et_tr_set = [];
       end
       
       %for vertical
       mod_data_name = 'full_eyetrack_initmods_vbars.mat';
       anal_name = 'full_eyetrack_vbars.mat';
       if exist(anal_name,'file') && use_ori(2)
           ver_data = load(anal_name,'it_mods','it_R2','dit_mods*','dit_R2*','et_tr_set');
           ver_mod_data = load(mod_data_name,'all_mod_SU*');
        else
           ver_mod_data.all_mod_SUnum = [];
           ver_data.et_tr_set = [];
       end
       
       clear cur_su_data
       hor_mu_numbers = hor_data.et_tr_set;
       ver_mu_numbers = ver_data.et_tr_set;
       hor_su_numbers = hor_mod_data.all_mod_SUnum(hor_mu_numbers);
       ver_su_numbers = ver_mod_data.all_mod_SUnum(ver_mu_numbers);
       if cur_expt_num ~= 86 | ~use_ori(1)
           hor_mu_numbers = hor_mu_numbers(hor_su_numbers > 0);
           ver_mu_numbers = ver_mu_numbers(ver_su_numbers > 0);
           hor_su_numbers = hor_su_numbers(hor_su_numbers > 0);
           ver_su_numbers = ver_su_numbers(ver_su_numbers > 0);
       end
       [tot_use_mus,IA] = unique([hor_mu_numbers ver_mu_numbers]);
       tot_use_sus = [hor_su_numbers ver_su_numbers];
       tot_use_sus = tot_use_sus(IA);
       for ii = 1:length(tot_use_mus)
           cur_hor_loc = hor_mu_numbers(hor_mu_numbers == tot_use_mus(ii));
           if isempty(cur_hor_loc);
               cur_hor_imp = -Inf;
           else
               cur_hor_imp = hor_data.dit_R2_LOO(ii,end,cur_hor_loc);
           end;
           cur_ver_loc = ver_mu_numbers(ver_mu_numbers == tot_use_mus(ii));
           if isempty(cur_ver_loc)
               cur_ver_imp = -Inf;
           else
               cur_ver_imp = ver_data.dit_R2_LOO(ii,end,cur_ver_loc);
           end
           if cur_hor_imp > cur_ver_imp
               cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 0;
               cur_su_data(ii).before_mod = hor_data.it_mods{1}(cur_hor_loc);
               cur_su_data(ii).after_mod = hor_data.dit_mods_LOO{ii,end}(cur_hor_loc);
               cur_su_data(ii).before_R2 = hor_data.it_R2(1,cur_hor_loc);
               cur_su_data(ii).after_R2 = cur_hor_imp;
               cur_su_data(ii).after_R2_fix = hor_data.it_R2(end,cur_hor_loc);
               cur_su_data(ii).probe_num = cur_probe_num; 
%                cur_su_data(ii).SU_num = 
           elseif cur_ver_imp > cur_hor_imp
               cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 90;
               cur_su_data(ii).before_mod = ver_data.it_mods{1}(cur_ver_loc);
               cur_su_data(ii).after_mod = ver_data.dit_mods_LOO{ii,end}(cur_ver_loc);
               cur_su_data(ii).before_R2 = ver_data.it_R2(1,cur_ver_loc);
               cur_su_data(ii).after_R2 = cur_ver_imp;
               cur_su_data(ii).after_R2_fix = ver_data.it_R2(end,cur_ver_loc);
               if cur_ver_loc > 96
                   cur_probe_num = ver_mod_data.all_mod_SU(cur_ver_loc);
                   cur_ver_spkloc = find(ver_SU_spk_numbers == tot_use_sus(ii));
                   cur_su_data(ii).num_blocks = cur_ver_SU_data.sua_data(cur_ver_spkloc).num_blocks;
                   cur_su_data(ii).n_spikes = cur_ver_SU_data.sua_data(cur_ver_spkloc).n_spikes;
                   cur_su_data(ii).mean_rate = cur_ver_SU_data.sua_data(cur_ver_spkloc).mean_rate;
               else
                   cur_probe_num = cur_ver_loc;
                   cur_su_data(ii).num_blocks = Inf;
                   cur_su_data(ii).n_spikes = cur_ver_SU_data.mua_data(cur_ver_loc).n_spikes;
                   cur_su_data(ii).mean_rate = cur_ver_SU_data.mua_data(cur_ver_loc).mean_rate;
               end
               cur_su_data(ii).probe_num = cur_probe_num;                             
           else
               cur_su_data(ii).usable = false;
               cur_su_data(ii).num_blocks = 0;
           end
           cur_su_data(ii).dx = 0.0565/2/scale_fac;
           cur_su_data(ii).monk_name = cur_expt_name(1);
           cur_su_data(ii).expt_num = cur_expt_num;
           cur_su_data(ii).rf_pos = rf_pos;
           cur_su_data(ii).rf_ecc = sqrt(sum(rf_pos.^2));
       end
       
       cur_su_data(~[cur_su_data(:).usable]) = [];
       
       fprintf('Using %d SUs\n',length(cur_su_data));
       all_su_data = cat(2,all_su_data,cur_su_data);
       
   else %FOR LEM EXPTS
%        anal_dir = ['~/Analysis/bruce/' cur_expt_name '/ET_final/'];
       anal_dir = ['~/Analysis/bruce/' cur_expt_name '/ET_final_imp/'];
       cd(anal_dir);
       
       clear dit_*
       mod_data_name = 'full_eyetrack_initmods.mat';
       anal_name = 'full_eyetrack_Cprior.mat';
       if exist(anal_name,'file')
           load(anal_name,'it_mods','it_R2','dit_mods*','dit_R2*','et_tr_set','et_params');
       end
           load(mod_data_name,'all_mod_SU*');
%        else
%            error('No data');
%        end
       clear cur_su_data
       mu_numbers = et_tr_set;
       su_numbers = all_mod_SUnum(mu_numbers);
       cur_SU_data = lem_SU_data{cellfun(@(x) x.Expt_num,lem_SU_data) == cur_expt_num};
       SU_spk_numbers = [cur_SU_data.sua_data(:).su_num];
      
       for ii = 1:length(mu_numbers)
           cur_loc = mu_numbers(ii);
           if ~isempty(cur_loc);
%                cur_imp = dit_R2(end,cur_loc);
               cur_imp = dit_R2_LOO(ii,end,cur_loc);
               cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 0;
               cur_su_data(ii).before_mod = it_mods{1}(cur_loc);
               cur_su_data(ii).after_mod = dit_mods_LOO{ii,end}(cur_loc);
               cur_su_data(ii).before_R2 = it_R2(1,cur_loc);
               cur_su_data(ii).after_R2 = cur_imp;
               cur_su_data(ii).after_R2_fix = it_R2(end,cur_loc);
               if cur_loc > 24
                   cur_probe_num = all_mod_SU(cur_loc);
                   cur_spkloc = find(SU_spk_numbers == su_numbers(ii));
                    cur_su_data(ii).num_blocks = cur_SU_data.sua_data(cur_spkloc).num_blocks;
                   cur_su_data(ii).n_spikes = cur_SU_data.sua_data(cur_spkloc).n_spikes;
                   cur_su_data(ii).mean_rate = cur_SU_data.sua_data(cur_spkloc).mean_rate;
              else
                  cur_probe_num = cur_loc;
                  cur_su_data(ii).num_blocks = Inf;
                  cur_su_data(ii).n_spikes = cur_SU_data.mua_data(cur_loc).n_spikes;
                  cur_su_data(ii).mean_rate = cur_SU_data.mua_data(cur_loc).mean_rate;
               end
               cur_su_data(ii).probe_num = cur_probe_num;
           else
               cur_su_data(ii).usable = false;
               cur_su_data(ii).num_blocks = 0;
           end
           %            if cur_expt_num < 240
           %                cur_su_data(ii).dx  = 0.125/4;
           % %            elseif cur_expt_num == 270 || cur_expt_num >= 280
           %            elseif cur_expt_num == 270 || cur_expt_num >= 280
           %                load(anal_name,'et_sp_dx');
           %                cur_su_data(ii).dx = et_sp_dx;
           %            else
           %                cur_su_data(ii).dx = 0.0565/2/scale_fac;
           %            end
           cur_su_data(ii).dx = et_params.sp_dx;
           cur_su_data(ii).monk_name = cur_expt_name(1);
           cur_su_data(ii).expt_num = cur_expt_num;
           cur_su_data(ii).rf_pos = rf_pos;
           cur_su_data(ii).rf_ecc = sqrt(sum(rf_pos.^2));
       end
       
       cur_su_data(~[cur_su_data(:).usable]) = [];
       
       fprintf('Using %d SUs\n',length(cur_su_data));
       all_su_data = cat(2,all_su_data,cur_su_data);
   end
    
end

%%
clear after_filt_data before_filt_data
for ii = 1:length(all_su_data)
    
        fprintf('Processing filters for unit %d of %d\n',ii,length(all_su_data));

%     rf_pos = all_su_data(ii).rf_pos;
%     bar_ori = all_su_data(ii).ori;
%     rf_proj = rf_pos(:,1)*cos((bar_ori-90)*pi/180) + rf_pos(:,2)*sin((bar_ori-90)*pi/180);
%     pix_ax = ((1:use_nPix_us)-use_nPix_us/2-0.5)*sp_dx + rf_proj;

    Xtargs = [all_su_data(ii).after_mod.mods(:).Xtarget];
    cor_filters = [all_su_data(ii).after_mod.mods(Xtargs==1).filtK];
    cor_stim_params = all_su_data(ii).after_mod.stim_params(1);
    after_filt_data(ii) = get_filter_properties(cor_filters,cor_stim_params,all_su_data(ii).dx);
 
%     Xtargs = [all_su_data(ii).before_mod.mods(:).Xtarget];
%     cor_filters = [all_su_data(ii).before_mod.mods(Xtargs==1).filtK];
%     cor_stim_params = all_su_data(ii).before_mod.stim_params(1);
%     before_filt_data(ii) = get_filter_properties(cor_filters,cor_stim_params,all_su_data(ii).dx);

end
%%
cd ~/Analysis/bruce/ET_final/
% save pooled_mu_data_v4test all_su_data after_filt_data
% save pooled_mu_data_newtemp all_su_data after_filt_data

%% COMPUTE SIMPLICITY
% n_frames = 1e5;
% flen = 12;
% simp = nan(length(all_su_data),1);
% weight_mean = nan(length(all_su_data),1);
% weight_std = nan(length(all_su_data),1);
% best_mean = nan(length(all_su_data),1);
% best_std = nan(length(all_su_data),1);
% 
% after_bsq_gests = reshape([after_filt_data(:).bsq_gest],[6 length(all_su_data)]);
% after_sq_gests = reshape([after_filt_data(:).sq_gest],[6 length(all_su_data)]);
% after_lin_gests = reshape([after_filt_data(:).lin_gest],[6 length(all_su_data)]);
% for ee = 1:n_expts
%     ee
%     su_set = find([all_su_data(:).expt_num] ==  str2num(expt_list{ee}(2:end)));
%     sp = all_su_data(su_set(1)).after_mod.stim_params(1);
%     dx = all_su_data(su_set(1)).dx;
%     cur_expt_num = str2num(expt_list{ee}(2:end));
%     if cur_expt_num > 281
%         spatial_usfac = 4;
%     else
%         spatial_usfac = 2;
%     end
%     
%     nPix = sp.stim_dims(2)/spatial_usfac;
%     if ismember(expt_list{ee},{'G091','G092','G093'})
%         dds = 0.67;
%     else
%         dds = 0.12;
%     end
%     rand_stim = zeros(n_frames,nPix);
%     rand_mat = rand(n_frames,nPix);
%     rand_stim(rand_mat <= dds/2) = -1;
%     rand_stim(rand_mat >= 1-(dds/2)) = 1;
%     
%     rand_stimmat_up = zeros(size(rand_stim,1),nPix*spatial_usfac);
%     for ii = 1:size(rand_stim,2)
%         for jj = 1:spatial_usfac
%             rand_stimmat_up(:,spatial_usfac*(ii-1)+jj) = rand_stim(:,ii);
%         end
%     end
%     nPix = nPix*spatial_usfac;
%     stim_params = NMMcreate_stim_params([flen nPix],0.01);
%     test_stim = create_time_embedding(rand_stimmat_up,stim_params);
%     
%     
%     for ii = 1:length(su_set)
%         Xtargs = [all_su_data(su_set(ii)).after_mod.mods(:).Xtarget];
%         cor_filters = [all_su_data(su_set(ii)).after_mod.mods(Xtargs==1).filtK];
%         filt_out = test_stim*cor_filters;
%         filt_out(:,2:end) = filt_out(:,2:end).^2;
%         filt_out_var = std(filt_out);
% %         filt_out_var = var(filt_out);
%         filt_out_var = filt_out_var/nansum(filt_out_var);
%         simp(su_set(ii)) = filt_out_var(1);
%         
% %         all_means = [after_filt_data(su_set(ii)).lin_mean after_filt_data(su_set(ii)).sq_mean ...
% %             after_filt_data(su_set(ii)).bsq_mean];
%         all_means = [after_lin_gests(1,su_set(ii)) after_sq_gests(1,su_set(ii)) ...
%             after_bsq_gests(1,su_set(ii))];
%         all_means = all_means - dx - nPix/2*dx; %center to deg relative to screen center
% 
%         all_sfs = [after_lin_gests(3,su_set(ii)) after_sq_gests(3,su_set(ii)) ...
%             after_bsq_gests(3,su_set(ii))];
% 
% %         all_stds = [after_filt_data(su_set(ii)).lin_std after_filt_data(su_set(ii)).sq_std ...
% %             after_filt_data(su_set(ii)).bsq_std];
%         all_stds = [after_lin_gests(2,su_set(ii)) after_sq_gests(2,su_set(ii)) ...
%             after_bsq_gests(2,su_set(ii))];
%         cur_ford = after_filt_data(su_set(ii)).ford;
%         weight_mean(su_set(ii)) = nansum(all_means.*filt_out_var(cur_ford));
%         weight_std(su_set(ii)) = nansum(all_stds.*filt_out_var(cur_ford));
%         weight_sf(su_set(ii)) = nansum(all_sfs.*filt_out_var(cur_ford));
%         [~,best_filt] = max(filt_out_var);
%         best_mean(su_set(ii)) = all_means(best_filt);
%         best_std(su_set(ii)) = all_stds(best_filt);
%     end
% end

n_frames = 1e5; %number of frames to create test stim
flen = 12;
simp = nan(length(all_su_data),1);
simp_before = nan(length(all_su_data),1);
weight_mean = nan(length(all_su_data),1);
weight_std = nan(length(all_su_data),1);
weight_sf = nan(length(all_su_data),1);

after_bsq_gests = reshape([after_filt_data(:).bsq_gest],[6 length(all_su_data)]);
after_sq_gests = reshape([after_filt_data(:).sq_gest],[6 length(all_su_data)]);
after_lin_gests = reshape([after_filt_data(:).lin_gest],[6 length(all_su_data)]);
for ee = 1:n_expts
    ee
    su_set = find([all_su_data(:).expt_num] ==  str2num(expt_list{ee}(2:end)));
    sp = all_su_data(su_set(1)).after_mod.stim_params(1);
    dx = all_su_data(su_set(1)).dx;
    cur_expt_num = str2num(expt_list{ee}(2:end));
    if cur_expt_num > 281
        spatial_usfac = 4;
    else
        spatial_usfac = 2;
    end

    %create a test stim with the appropriate stats to evaluate relative
    %filter outputs
     nPix = sp.stim_dims(2)/spatial_usfac; %number of bars
%    nPix = sp.stim_dims(2);
    if ismember(expt_list{ee},{'G091','G092','G093'})
        dds = 0.67;
    else
        dds = 0.12;
    end
    
    %create random bar stim
    rand_stim = zeros(n_frames,nPix);
    rand_mat = rand(n_frames,nPix);
    rand_stim(rand_mat <= dds/2) = -1;
    rand_stim(rand_mat >= 1-(dds/2)) = 1;
%     stim_params = NMMcreate_stim_params([flen nPix],0.01);
%     test_stim = create_time_embedding(rand_stim,stim_params);    

%up-sample bars
        rand_stimmat_up = zeros(size(rand_stim,1),nPix*spatial_usfac);
    for ii = 1:size(rand_stim,2)
        for jj = 1:spatial_usfac
            rand_stimmat_up(:,spatial_usfac*(ii-1)+jj) = rand_stim(:,ii);
        end
    end
    nPix = nPix*spatial_usfac;
    stim_params = NMMcreate_stim_params([flen nPix],0.01);
    test_stim = create_time_embedding(rand_stimmat_up,stim_params);

    for ii = 1:length(su_set)
        Xtargs = [all_su_data(su_set(ii)).after_mod.mods(:).Xtarget];
                
%repeat for corrected model
        cor_filters = [all_su_data(su_set(ii)).after_mod.mods(Xtargs==1).filtK];
        offset = all_su_data(su_set(ii)).after_mod.spk_NL_params(1);
        beta = all_su_data(su_set(ii)).after_mod.spk_NL_params(2);
        alpha = all_su_data(su_set(ii)).after_mod.spk_NL_params(3);
        offset = offset + mean(all_su_data(su_set(ii)).after_mod.mods(4).filtK);
        filt_out = test_stim*cor_filters;
        filt_out(:,2:end) = filt_out(:,2:end).^2;
         rate_out = alpha*log(1+exp(beta*(sum(filt_out,2) + offset)));
         
         filt_out_var = std(filt_out);
         filt_out_var = filt_out_var/nansum(filt_out_var);
         simp(su_set(ii)) = filt_out_var(1);
         
        filt_outr = -test_stim*cor_filters;
        filt_outr(:,2:end) = filt_outr(:,2:end).^2;
        rate_out_rev = alpha*log(1+exp(beta*(sum(filt_outr,2) + offset)));
        
        %phase-reversal modulation for corrected model
        prm_after(su_set(ii)) = mean(abs(rate_out_rev - rate_out))/mean(rate_out);

        all_means = [after_lin_gests(1,su_set(ii)) after_sq_gests(1,su_set(ii)) ...
            after_bsq_gests(1,su_set(ii))];
        all_means = all_means - dx - nPix/2*dx; %center to deg relative to screen center

        all_sfs = [after_lin_gests(3,su_set(ii)) after_sq_gests(3,su_set(ii)) ...
            after_bsq_gests(3,su_set(ii))];
        
        all_stds = [after_lin_gests(2,su_set(ii)) after_sq_gests(2,su_set(ii)) ...
            after_bsq_gests(2,su_set(ii))];
        cur_ford = after_filt_data(su_set(ii)).ford;
        weight_mean(su_set(ii)) = nansum(all_means.*filt_out_var(cur_ford));
        weight_std(su_set(ii)) = nansum(all_stds.*filt_out_var(cur_ford));
        weight_sf(su_set(ii)) = nansum(all_sfs.*filt_out_var(cur_ford));
    end
end

%%
num_blocks = [all_su_data(:).num_blocks];
expt_nums = [all_su_data(:).expt_num];
probe_nums = [all_su_data(:).probe_num];
n_spikes = [all_su_data(:).n_spikes];
mean_rate = [all_su_data(:).mean_rate];
expt_oris = [all_su_data(:).ori];
before_R2 = [all_su_data(:).before_R2];
after_R2 = [all_su_data(:).after_R2];

%% compute RF positions
load ~/Data/bruce/general_array_data/array_pos_data.mat
interp_ecc = sqrt(interp_x.^2 + interp_y.^2);

all_ecc = interp_ecc(probe_nums);
all_Y = interp_y(probe_nums);
all_X = interp_x(probe_nums);

lemExpts = [266 270 275 277 281 287 289 294];
lemUnits = ismember(expt_nums,lemExpts);
jbeUnits = ~ismember(expt_nums,lemExpts);
expt_oris(expt_nums == 266) = 80;
expt_oris(expt_nums == 270) = 60;
expt_oris(expt_nums == 275) = 135;
expt_oris(expt_nums == 277) = 70;
expt_oris(expt_nums == 281) = 140;
expt_oris(expt_nums == 287) = 90;
expt_oris(expt_nums == 289) = 160;
expt_oris(expt_nums == 294) = 40;

clear rf_par rf_orth
rf_pos = reshape([all_su_data(:).rf_pos],[2 length(all_su_data)]);


rf_orth = [weight_mean']; %component of RF position orthoganol to bar stim
rf_orth_avg = -rf_pos(1,:).*sind(expt_oris) - rf_pos(2,:).*cosd(expt_oris); %avg RF position (in eye coords, so flip hori comp for numbers stored in screen coords)
rf_orth = rf_orth + rf_orth_avg;

%for RF component parallel to bar stim, take avg RF position estimates for
%each probe
rf_par(jbeUnits) = -all_X(jbeUnits).*cosd(expt_oris(jbeUnits)') + all_Y(jbeUnits).*sind(expt_oris(jbeUnits)');
rf_par(lemUnits) = -rf_pos(1,lemUnits)'.*cosd(expt_oris(lemUnits)') + rf_pos(2,lemUnits)'.*sind(expt_oris(lemUnits)');

rf_eccs = sqrt(rf_par.^2 + rf_orth.^2);

%%
all_unit_data.num_blocks = num_blocks;
all_unit_data.expt_nums = expt_nums;
all_unit_data.probe_nums = probe_nums;
all_unit_data.n_spikes = n_spikes;
all_unit_data.mean_rate = mean_rate;
all_unit_data.expt_oris = expt_oris;
all_unit_data.before_R2 = before_R2;
all_unit_data.after_R2 = after_R2;
all_unit_data.rf_eccs = rf_eccs;
all_unit_data.rf_pos = rf_pos;
all_unit_data.prm = prm_after;
all_unit_data.rf_width = weight_std*2;
all_unit_data.rf_SF = weight_sf;

cd /home/james/Analysis/bruce/ET_final
save unit_summary_data2 all_unit_data
%%
min_num_blocks = 3;
use_R2_thresh = 0.00; %minimum post-correction R2 value to include the unit

% used_units = find(num_blocks >= min_num_blocks);
used_units = find(num_blocks >= min_num_blocks & after_R2 >= use_R2_thresh);
badu = [];
% badu = [badu find([all_su_data(:).expt_num] == 266 & [all_su_data(:).probe_num] == 24)]; %this deepest probe is either in V2 or a different patch of V1
badu = [badu find(expt_nums == 270 & [all_su_data(:).probe_num] == 22,1)]; %this MU is picking up on two spatially distinct RFs and so the Gabor estimation fails miserably.
used_units(ismember(used_units,badu)) = [];

after_R2 = after_R2(used_units);
before_R2 = before_R2(used_units);
R2_imp_fac = after_R2./before_R2;
RF_widths = weight_std(used_units)*2;
rf_eccs = rf_eccs(used_units);
expt_nums = expt_nums(used_units);
num_blocks = num_blocks(used_units);

after_sq_gests = reshape([after_filt_data(used_units).sq_gest],[6 length(used_units)]);
after_lin_gests = reshape([after_filt_data(used_units).lin_gest],[6 length(used_units)]);

thresh_R2 = 0.8;
after_sq_gR2 = [after_filt_data(used_units).sq_est_R2];
after_lin_gR2 = [after_filt_data(used_units).lin_est_R2];

%%
close all
check_expts = [270];
cur_set = find(ismember(expt_nums,check_expts));
for ii = 1:length(cur_set)
   NMMdisplay_model(all_su_data(used_units(cur_set(ii))).after_mod,[],[],1)
   fprintf('Expt %d, Probe: %d\n',expt_nums(cur_set(ii)),all_su_data(used_units(cur_set(ii))).probe_num);
   fprintf('R2 imp: %.4f\n\n',R2_imp_fac(cur_set(ii)));
   fprintf('Ov width %.4f, lin: %.4f Sq: %.4f\n\n',weight_std(used_units(cur_set(ii)))*2,after_lin_gests(2,cur_set(ii))*2,after_sq_gests(2,cur_set(ii))*2);
   
%    figure
%    plot(after_filt_data(used_units(cur_set(ii))).sq_gabor_est); hold on 
%    plot(after_filt_data(used_units(cur_set(ii))).sq_slice,'r'); hold on 
   pause
   close
end

%% RF width vs R2 improvement WEIGHTED WIDTH
% rmpath('~/James_scripts/bruce/bruce_code/')
 
fig_width = 3.27; rel_height = 0.8;
mSize = 10;

% allUnits = find(after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2);
allUnits = find((after_sq_gR2 >= thresh_R2 | after_lin_gR2 >= thresh_R2) & ~isinf(num_blocks));
% allUnits = 1:length(used_units);

% lemExpts = [266 270 275 277 281 287 289 294];
lemExpts = [266 275 ];
lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
lemExpts2 = [270 277 281 287 289 294];
lemUnits2 = allUnits(ismember(expt_nums(allUnits),lemExpts2));
lemExpts = [266 270 275 277 281 287 289 294];
jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));

X = RF_widths;
Y = R2_imp_fac;

% h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth.fig');
h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth3.fig');
hold on
plot(X(jbeUnits),Y(jbeUnits),'.','markersize',mSize);
plot(X(lemUnits),Y(lemUnits),'r.','markersize',mSize);
plot(X(lemUnits2),Y(lemUnits2),'g.','markersize',mSize);
line_hand = findall(gcf,'type','line','-and','color','k');
uistack(line_hand,'top');
figufy(h);

% xlim([0.055 0.6]);
xlim([0.08 0.85]);
% ylim([1 5]);
ylim([1 6]);
set(gca,'xscale','log','yscale','log');

% F = @(x,xdata)x(1)*xdata.^(x(2));
% x0 = [4 -0.5];
% [x,resnorm] = lsqcurvefit(F,x0,X,Y');
% xx = linspace(0.01,1,100);
% plot(xx,x(1)*xx.^x(2),'g','linewidth',2);
% 
% fname = [fig_dir 'modimp_vs_RFwidth_overlay9.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);
% 

%% Simplicity vs R2 improvement 
% fig_width = 3.27; rel_height = 0.8;
% mSize = 10;
% 
% % allUnits = find(after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2);
% allUnits = find(after_sq_gR2 >= thresh_R2 | after_lin_gR2 >= thresh_R2);
% % allUnits = 1:length(used_units);
% 
% lemExpts = [266 270 275 277 281 287 289 294];
% lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
% jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));
% 
% X = simp(used_units);
% Y = R2_imp_fac;
% 
% h1 = figure(); hold on
% plot(X(jbeUnits),Y(jbeUnits),'.','markersize',mSize*0.6);
% plot(X(lemUnits),Y(lemUnits),'r.','markersize',mSize);
% 
% mdl = LinearModel.fit(X,Y,'linear','Robustopts','on');
% xr = linspace(nanmin(X),nanmax(X),100);
% [ypred,yci] = predict(mdl,xr(:));
% plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
% plot(xr,yci,'--','color',[0.2 0.8 0.2]);
% yl = ylim();
% ylim([1 6.1]);
% set(gca,'yscale','log')
% 
% fname = [fig_dir 'modimp_vs_simp.pdf'];
% % exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h1);

%% RF width vs eccentricity
mSize = 10;

% allUnits = find(after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2);
allUnits = find(after_sq_gR2 >= thresh_R2 | after_lin_gR2 >= thresh_R2);
% allUnits = 1:length(used_units);

% lemExpts = [266 270 275 277 281 287 289 294];
lemExpts = [266 275];
lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
lemExpts2 = [270 277 281 287 289 294];
lemUnits2 = allUnits(ismember(expt_nums(allUnits),lemExpts2));
lemExpts = [266 270 275 277 281 287 289 294];
jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));

jit_amp = 0;
X = rf_eccs + randn(size(rf_eccs))*jit_amp;
Y = RF_widths;

h = figure(); hold on
plot(X(jbeUnits),Y(jbeUnits),'.','markersize',mSize);
plot(X(lemUnits),Y(lemUnits),'r.','markersize',mSize);
plot(X(lemUnits2),Y(lemUnits2),'g.','markersize',mSize);

mdl = LinearModel.fit(X(allUnits),Y(allUnits),'linear','Robustopts','on');
% mdl = LinearModel.fit(X(allUnits),Y(allUnits),'linear','Robustopts','off');
xr = linspace(nanmin(X(allUnits)),nanmax(X(allUnits)),100);
[ypred,yci] = predict(mdl,xr(:));
plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
plot(xr,yci,'--','color',[0.2 0.8 0.2]);
xlabel('RF ecc (deg)');
ylabel('RF width (deg)');
set(gca,'yscale','log');
% ylim([0.05 0.5]);
ylim([0.08 0.8]);

fig_width = 3.27; rel_height = 0.8;
figufy(h);
% fname = [fig_dir 'RFwidth_vs_ecc8.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% SF vs eccentricity
lemExpts = [266 270 275 277];

lemUnits = ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
jbeUnits = ~ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
allUnits = find(after_lin_gR2 >= thresh_R2);

% jit_amp = 0.01;
jit_amp = 0.0;
X = rf_eccs + randn(size(rf_eccs))*jit_amp;
Y = after_lin_gests(3,:);

h = figure(); hold on
plot(X(jbeUnits),Y(jbeUnits),'.');
plot(X(lemUnits),Y(lemUnits),'r.');

mdl = LinearModel.fit(X(allUnits),Y(allUnits),'linear','Robustopts','on');
xr = linspace(nanmin(X(allUnits)),nanmax(X(allUnits)),100);
[ypred,yci] = predict(mdl,xr(:));
plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
plot(xr,yci,'--','color',[0.2 0.8 0.2]);
xlabel('RF ecc (deg)');
ylabel('SF (1/deg)');
% set(gca,'yscale','log');
% ylim([0.05 0.5]);

fig_width = 3.27; rel_height = 0.8;
figufy(h);
% fname = [fig_dir 'RFwidth_vs_ecc.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% Width vs SF
lemExpts = [266 270 275 277];

lemUnits = ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
jbeUnits = ~ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
allUnits = find(after_lin_gR2 >= thresh_R2);

X = RF_widths;
Y = after_lin_gests(3,:);

h = figure(); hold on
plot(X(jbeUnits),Y(jbeUnits),'.');
plot(X(lemUnits),Y(lemUnits),'r.');

mdl = LinearModel.fit(X(allUnits),Y(allUnits),'linear','Robustopts','on');
xr = linspace(nanmin(X(allUnits)),nanmax(X(allUnits)),100);
[ypred,yci] = predict(mdl,xr(:));
plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
plot(xr,yci,'--','color',[0.2 0.8 0.2]);
xlabel('RF width (deg)');
ylabel('SF (1/deg)');
% set(gca,'yscale','log');
% ylim([0.05 0.5]);

fig_width = 3.27; rel_height = 0.8;
figufy(h);
fname = [fig_dir 'RFwidth_vs_ecc2.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);


%%
% all_RF_widths = after_sq_gests(2,:)*2;
% all_RF_widths(after_sq_gR2 < thresh_R2) = nan;
% save('~/Analysis/bruce/ET_final/G086_hbar_mua_RFwidth.mat','all_RF_widths');
%%
% X = rf_eccs;
% X = [after_filt_data(:).lin_FFx];
X = 2*after_sq_gests(2,:);
% X = [after_filt_data(:).simplicity];
% X = [after_filt_data(:).lin_std];
% Y = [after_filt_data(:).sq_std];
% Y = after_lin_gests(2,:);
Y = R2_imp_fac;

bad_pts = after_sq_gR2 < thresh_R2;
% bad_pts = after_lin_gR2 < thresh_R2;
% bad_pts = after_sq_gR2 < thresh_R2 | after_lin_gR2 < thresh_R2;
X(bad_pts) = nan;
Y(bad_pts) = nan;



figure;
clear g_*
hold on
for ii = 1:length(unique_expt_nums)
    g_avgX(ii) = nanmean(X(G==ii));
    g_semX(ii) = nanstd(X(G==ii))/sqrt(nansum(G==ii));
    g_avgY(ii) = nanmean(Y(G==ii));
    g_semY(ii) = nanstd(Y(G==ii))/sqrt(nansum(G==ii));
    cur_enum(ii) = unique(expt_nums(G==ii));
    if ismember(cur_enum(ii),sbar_set)
%         errorbarxy(g_avgX(ii),g_avgY(ii),g_semX(ii),g_semY(ii),{'r','r','r'});
        %         g_avgX(ii) = nan;
        %         X(G==ii) = nan;
        g_isS(G==ii) = true;
    else
        g_isS(G==ii) = false;
%         errorbarxy(g_avgX(ii),g_avgY(ii),g_semX(ii),g_semY(ii),{'k','k','k'});
    end
end
plot(X(~g_isS),Y(~g_isS),'b.');
% plot(g_avgX(~ismember(cur_enum,sbar_set)),g_avgY(~ismember(cur_enum,sbar_set)),'ko','linewidth',2,'markersize',8);
% plot(g_avgX(ismember(cur_enum,sbar_set)),g_avgY(ismember(cur_enum,sbar_set)),'go','linewidth',2);
plot(X(g_isS),Y(g_isS),'r.');
mdl = LinearModel.fit(X,Y,'linear','Robustopts','on')
xr = linspace(nanmin(X),nanmax(X),100);
% [ypred,yci] = predict(mdl,xr(:));
% plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
% plot(xr,yci,'--','color',[0.2 0.8 0.2]);
% aoctool(X,Y,G,.05,[],[],[],'on','parallel lines');

%%
% plot([before_filt_data(:).lin_FFx],[after_filt_data(:).lin_FFx],'o');
% plot([before_filt_data(:).lin_tot_pow],[after_filt_data(:).lin_tot_pow],'o');
% 
% plot([after_filt_data(:).sq_std],R2_imp_fac,'o');
% 
% plot([all_su_data(:).rf_ecc],R2_imp_fac,'o');
% plot([all_su_data(:).rf_ecc],[after_filt_data(:).sq_FFx],'o');
% 

%% RF width vs R2 improvement WEIGHTED WIDTH
% rmpath('~/James_scripts/bruce/bruce_code/')
 
fig_width = 3.27; rel_height = 0.8;
mSize = 10;

allUnits = find(after_sq_gR2 >= 0.5);

lemExpts = [266 270 275 277 281 287 289 294];
lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));

% h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth.fig');
h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth3.fig');
hold on
plot(after_sq_gests(2,jbeUnits)*2,R2_imp_fac(jbeUnits),'.','markersize',mSize);
plot(after_sq_gests(2,lemUnits)*2,R2_imp_fac(lemUnits),'r.','markersize',mSize);
line_hand = findall(gcf,'type','line','-and','color','k');
uistack(line_hand,'top');
figufy(h);

% X = 2*weight_std(used_units);
% % X = 2*after_sq_gests(2,:);
% Y = R2_imp_fac;
% % Y = R2_imp_fac_fix;
% n_bins = 20;
% pbin_edges = prctile(X,linspace(0,100,n_bins+1));
% avg_X = nan(n_bins,1);
% avg_Y = nan(n_bins,1);
% std_X = nan(n_bins,1);
% std_Y = nan(n_bins,1);
% for ii = 1:n_bins
%    curset = find(X >= pbin_edges(ii) & X < pbin_edges(ii+1));
%    avg_X(ii) = mean(X(curset));
%    avg_Y(ii) = mean(Y(curset));
%    std_X(ii) = std(X(curset))/sqrt(length(curset));
%    std_Y(ii) = std(Y(curset))/sqrt(length(curset));
% %    std_X(ii) = std(X(curset));
% %    std_Y(ii) = std(Y(curset));
% end

% eh = errorbarxy(avg_X, avg_Y, std_X, std_Y,{'g' 'g' 'g'});
% set(eh(1),'linewidth',2);
% set(eh(2:end),'linewidth',2);

% [~,xord] = sort(X(allUnits));
% ysm = smooth(Y(allUnits(xord)),50,'rloess');
% plot(X(allUnits(xord)),ysm,'g-','linewidth',2);

% xlim([0.055 0.6]);
xlim([0.07 0.8]);
% ylim([1 5]);
ylim([1 6]);
set(gca,'xscale','log','yscale','log');

% F = @(x,xdata)x(1)*xdata.^(x(2));
% x0 = [4 -0.5];
% [x,resnorm] = lsqcurvefit(F,x0,X,Y');
% xx = linspace(0.01,1,100);
% plot(xx,x(1)*xx.^x(2),'g','linewidth',2);
% 
% fname = [fig_dir 'modimp_vs_RFwidth_overlay6.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);
