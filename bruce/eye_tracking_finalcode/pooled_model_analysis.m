clear all
close all

addpath('~/James_scripts/bruce/eye_tracking_finalcode/');
% expt_list = {'G086'};
% expt_list = {'M270'};
expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
% expt_list = {'M232','M235','M239'};
% expt_list = {'M266','M270','M275','M277'};
% expt_list = {'M266','M270','M275','M277','G085','G086','G087','G088','G089','G091','G093','G095'};
% expt_list = {'M232','M235','M239','M266','M270','M275','M277','G085','G086','G087','G088','G089','G091','G093','G095'};

print_out = false;

fig_dir = '/home/james/Analysis/bruce/ET_final/';

% use_ori = [true true];
use_ori = [true false];
% use_ori = [false true];

hor_SU_data = load('~/Analysis/bruce/ET_final/SUA_sum_data_hor.mat');
hor_SU_data = hor_SU_data.Expt_sua_data;

ver_SU_data = load('~/Analysis/bruce/ET_final/SUA_sum_data_ver.mat');
ver_SU_data = ver_SU_data.Expt_sua_data;

lem_SU_data = load('~/Analysis/bruce/ET_final/SUA_sum_data_lem.mat');
lem_SU_data = lem_SU_data.Expt_sua_data;

%%
n_expts = length(expt_list);
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
    
    data_dir = ['~/Data/bruce/' cur_expt_name];
    if cur_expt_name(1) == 'G'
        load([data_dir sprintf('/jbe%sExpts.mat',cur_expt_name)]); %load in Expts struct
    else
        load([data_dir sprintf('/lem%sExpts.mat',cur_expt_name)]); %load in Expts struct
    end
    rf_pos = Expts{1}.Stimvals.rf(1:2)/scale_fac;
       
   if cur_expt_name(1) == 'G'       
       anal_dir = ['~/Analysis/bruce/' cur_expt_name '/ET_final/'];
       cd(anal_dir);
       
        clear dit_*
      %for horizontal
       mod_data_name = 'monoc_eyecorr_hbar_mods.mat';
%        anal_name = 'monoc_eyecorr_hbar.mat';
       anal_name = 'monoc_eyecorr_hbar2.mat';
       if exist(anal_name,'file') && use_ori(1)
           hor_data = load(anal_name,'it_mods','it_R2','dit_mods_LOO','dit_R2_LOO','it_R2_LOO','dit_R2','et_tr_set');
           hor_mod_data = load(mod_data_name,'all_mod_SU*');
           cur_hor_SU_data = hor_SU_data{cellfun(@(x) x.Expt_num,hor_SU_data) == cur_expt_num};
        hor_SU_spk_numbers = [cur_hor_SU_data.sua_data(:).su_num];
      else
           hor_mod_data.all_mod_SUnum = [];
           hor_data.et_tr_set = [];
           cur_hor_SU_data = [];
           hor_SU_spk_numbers = [];
       end
       
       %for vertical
       mod_data_name = 'monoc_eyecorr_vbar_mods.mat';
%        anal_name = 'monoc_eyecorr_vbar.mat';
       anal_name = 'monoc_eyecorr_vbar2.mat';
       if exist(anal_name,'file') && use_ori(2)
           ver_data = load(anal_name,'it_mods','it_R2','dit_mods_LOO','dit_R2_LOO','it_R2_LOO','dit_R2','et_tr_set');
           ver_mod_data = load(mod_data_name,'all_mod_SU*');
           cur_ver_SU_data = ver_SU_data{cellfun(@(x) x.Expt_num,ver_SU_data) == cur_expt_num};
       ver_SU_spk_numbers = [cur_ver_SU_data.sua_data(:).su_num];
       else
           ver_mod_data.all_mod_SUnum = [];
           ver_data.et_tr_set = [];
           cur_ver_SU_data = [];
           ver_SU_spk_numbers = [];
       end
       
       clear cur_su_data
       hor_su_numbers = hor_mod_data.all_mod_SUnum(hor_data.et_tr_set(hor_mod_data.all_mod_SUnum(hor_data.et_tr_set) > 0));
       ver_su_numbers = ver_mod_data.all_mod_SUnum(ver_data.et_tr_set(ver_mod_data.all_mod_SUnum(ver_data.et_tr_set) > 0));
       tot_use_sus = unique([hor_su_numbers ver_su_numbers]);
       hor_su_inds = find(ismember(hor_mod_data.all_mod_SUnum,hor_su_numbers));
       ver_su_inds = find(ismember(ver_mod_data.all_mod_SUnum,ver_su_numbers));
       for ii = 1:length(tot_use_sus)
           cur_hor_loc = find(hor_su_numbers == tot_use_sus(ii));
           cur_hor_spkloc = find(hor_SU_spk_numbers == tot_use_sus(ii));
           if isempty(cur_hor_loc);
               cur_hor_imp = -Inf;
           else
               cur_hor_imp = hor_data.dit_R2_LOO(cur_hor_loc,end,hor_su_inds(cur_hor_loc));
               cur_hor_imp_fix = hor_data.it_R2_LOO(cur_hor_loc,end,hor_su_inds(cur_hor_loc));
           end;
           
           cur_ver_loc = find(ver_su_numbers == tot_use_sus(ii));
           cur_ver_spkloc = find(ver_SU_spk_numbers == tot_use_sus(ii));
           if isempty(cur_ver_loc)
               cur_ver_imp = -Inf;
           else
               cur_ver_imp = ver_data.dit_R2_LOO(cur_ver_loc,end,ver_su_inds(cur_ver_loc));
               cur_ver_imp_fix = ver_data.it_R2_LOO(cur_ver_loc,end,ver_su_inds(cur_ver_loc));
           end
           if cur_hor_imp > cur_ver_imp
            cur_hor_R2fract = hor_data.it_R2(1,hor_su_inds(cur_hor_loc))/sum(hor_data.it_R2(1,:));
              cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 0;
               cur_su_data(ii).before_mod = hor_data.it_mods{1}(hor_su_inds(cur_hor_loc));
               cur_su_data(ii).after_mod = hor_data.dit_mods_LOO{cur_hor_loc,end}(hor_su_inds(cur_hor_loc));
               cur_su_data(ii).R2_frac = cur_hor_R2fract;
               cur_su_data(ii).before_R2 = hor_data.it_R2(1,hor_su_inds(cur_hor_loc));
               cur_su_data(ii).after_R2 = cur_hor_imp;
               cur_su_data(ii).after_R2_fix = cur_hor_imp_fix;
               cur_su_data(ii).after_R2_fix_noXV = hor_data.it_R2(end,hor_su_inds(cur_hor_loc));
               cur_su_data(ii).after_R2_noXV = hor_data.dit_R2(end,hor_su_inds(cur_hor_loc));
               cur_su_data(ii).probe_num = hor_mod_data.all_mod_SU(hor_su_inds(cur_hor_loc));   
               cur_su_data(ii).num_blocks = cur_hor_SU_data.sua_data(cur_hor_spkloc).num_blocks;
               cur_su_data(ii).n_spikes = cur_hor_SU_data.sua_data(cur_hor_spkloc).n_spikes;
               cur_su_data(ii).mean_rate = cur_hor_SU_data.sua_data(cur_hor_spkloc).mean_rate;
           elseif cur_ver_imp > cur_hor_imp
             cur_ver_R2fract = ver_data.it_R2(1,ver_su_inds(cur_ver_loc))/sum(ver_data.it_R2(1,:));
             cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 90;
               cur_su_data(ii).before_mod = ver_data.it_mods{1}(ver_su_inds(cur_ver_loc));
               cur_su_data(ii).after_mod = ver_data.dit_mods_LOO{cur_ver_loc,end}(ver_su_inds(cur_ver_loc));
               cur_su_data(ii).R2_frac = cur_ver_R2fract;
               cur_su_data(ii).before_R2 = ver_data.it_R2(1,ver_su_inds(cur_ver_loc));
               cur_su_data(ii).after_R2 = cur_ver_imp;
               cur_su_data(ii).after_R2_fix = cur_ver_imp_fix;
               cur_su_data(ii).after_R2_fix_noXV = ver_data.it_R2(end,ver_su_inds(cur_ver_loc));
               cur_su_data(ii).after_R2_noXV = ver_data.dit_R2(end,ver_su_inds(cur_ver_loc));
               cur_su_data(ii).probe_num = ver_mod_data.all_mod_SU(ver_su_inds(cur_ver_loc));                             
               cur_su_data(ii).num_blocks = cur_ver_SU_data.sua_data(cur_ver_spkloc).num_blocks;
               cur_su_data(ii).n_spikes = cur_ver_SU_data.sua_data(cur_ver_spkloc).n_spikes;
               cur_su_data(ii).mean_rate = cur_ver_SU_data.sua_data(cur_ver_spkloc).mean_rate;
           else
               cur_su_data(ii).usable = false;
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
       anal_dir = ['~/Analysis/bruce/' cur_expt_name '/ET_final/'];
       cd(anal_dir);
       
       clear dit_*
       mod_data_name = 'monoc_eyecorr_mods.mat';
%        anal_name = 'monoc_eyecorr_Cprior.mat';
       anal_name = 'monoc_eyecorr2_Cprior.mat';
%         anal_name = 'monoc_eyecorr.mat';
      if exist(anal_name,'file')
           load(anal_name,'it_mods','it_R2','dit_mods_LOO','dit_R2_LOO','it_R2_LOO','dit_R2','et_tr_set');
           load(mod_data_name,'all_mod_SU*');
       else
           error('No data');
       end
       clear cur_su_data
       su_numbers = all_mod_SUnum(et_tr_set(all_mod_SUnum(et_tr_set) > 0));
       su_inds = find(ismember(all_mod_SUnum,su_numbers));
       
       cur_SU_data = lem_SU_data{cellfun(@(x) x.Expt_num,lem_SU_data) == cur_expt_num};
       SU_spk_numbers = [cur_SU_data.sua_data(:).su_num];
       for ii = 1:length(su_numbers)
           cur_loc = find(su_numbers == su_numbers(ii));
           cur_spkloc = find(SU_spk_numbers == su_numbers(ii));
           cur_R2fract = it_R2(1,su_inds(cur_loc))/sum(it_R2(1,:));
          if ~isempty(cur_loc);
               cur_imp = dit_R2_LOO(cur_loc,end,su_inds(cur_loc));
               cur_imp_fix = it_R2_LOO(cur_loc,end,su_inds(cur_loc));
               cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 0;
               cur_su_data(ii).before_mod = it_mods{1}(su_inds(cur_loc));
               cur_su_data(ii).after_mod = dit_mods_LOO{cur_loc,end}(su_inds(cur_loc));
               cur_su_data(ii).R2_frac = cur_R2fract;
               cur_su_data(ii).before_R2 = it_R2(1,su_inds(cur_loc));
               cur_su_data(ii).after_R2 = cur_imp;
               cur_su_data(ii).after_R2_fix = cur_imp_fix;
               cur_su_data(ii).after_R2_fix_noXV = it_R2(end,su_inds(cur_loc));
               cur_su_data(ii).after_R2_noXV = dit_R2(end,su_inds(cur_loc));
               cur_su_data(ii).probe_num = all_mod_SU(su_inds(cur_loc));
                cur_su_data(ii).num_blocks = cur_SU_data.sua_data(cur_spkloc).num_blocks;
               cur_su_data(ii).n_spikes = cur_SU_data.sua_data(cur_spkloc).n_spikes;
               cur_su_data(ii).mean_rate = cur_SU_data.sua_data(cur_spkloc).mean_rate;
          else
               cur_su_data(ii).usable = false;
           end
           if cur_expt_num < 240
               cur_su_data(ii).dx  = 0.125/4;
           else
               cur_su_data(ii).dx = 0.0565/2/scale_fac;
           end
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
min_num_blocks = 3;
num_blocks = [all_su_data(:).num_blocks];
n_spikes = [all_su_data(:).n_spikes];

used_sus = find(num_blocks >= min_num_blocks);

expt_nums = [all_su_data(used_sus).expt_num];
rf_eccs = [all_su_data(used_sus).rf_ecc];
lemExpts = [266 270 275 277];
lemUnits = ismember(expt_nums,lemExpts);
jbeUnits = ~ismember(expt_nums,lemExpts);


R2_fract = [all_su_data(used_sus).R2_frac];
after_R2 = [all_su_data(used_sus).after_R2];
after_R2_fix = [all_su_data(used_sus).after_R2_fix];
after_R2_fix_noXV = [all_su_data(used_sus).after_R2_fix_noXV];
before_R2 = [all_su_data(used_sus).before_R2];
after_R2_noXV = [all_su_data(used_sus).after_R2_noXV];
R2_imp_fac = after_R2./before_R2;
R2_imp_fac_noXV = after_R2_noXV./before_R2;
R2_imp_fac_fix = after_R2_fix./before_R2;
R2_imp_fac_fix_noXV = after_R2_fix_noXV./before_R2;

percent_diff_LOO = (R2_imp_fac_noXV - R2_imp_fac)./R2_imp_fac;
percent_diff_LOO_fix = (R2_imp_fac_fix_noXV - R2_imp_fac_fix)./R2_imp_fac_fix;

%% PLOT XV VS NON XV R2
mSize = 10;

h = figure; hold on
% plot(R2_imp_fac_noXV,R2_imp_fac,'.','markersize',mSize);
plot(after_R2_noXV,after_R2,'k.','markersize',mSize);
xlim([0 0.3]); ylim([0 0.3]);
xlabel('R2');
ylabel('Xval R2');
line([0 0.3],[0 0.3],'color','r');

fig_width = 3; 
rel_height = 0.85;

% figufy(h);
% fname = [fig_dir 'R2_vs_XVR2.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);


%%
close all
clear after_filt_data before_filt_data
for ii = 1:length(all_su_data)
    fprintf('Processing filters for unit %d of %d\n',ii,length(all_su_data));
    
    Xtargs = [all_su_data(ii).after_mod.mods(:).Xtarget];
    cor_filters = [all_su_data(ii).after_mod.mods(Xtargs==1).filtK];
    cor_stim_params = all_su_data(ii).after_mod.stim_params(1);
    after_filt_data(ii) = get_filter_properties(cor_filters,cor_stim_params,all_su_data(ii).dx);
 
    Xtargs = [all_su_data(ii).before_mod.mods(:).Xtarget];
    cor_filters = [all_su_data(ii).before_mod.mods(Xtargs==1).filtK];
    cor_stim_params = all_su_data(ii).before_mod.stim_params(1);
    before_filt_data(ii) = get_filter_properties(cor_filters,cor_stim_params,all_su_data(ii).dx);

end

%% COMPUTE SIMPLICITY
n_frames = 1e5;
flen = 12;
simp = nan(length(all_su_data),1);
simp_before = nan(length(all_su_data),1);
weight_mean = nan(length(all_su_data),1);
weight_std = nan(length(all_su_data),1);

after_bsq_gests = reshape([after_filt_data(:).bsq_gest],[6 length(all_su_data)]);
after_sq_gests = reshape([after_filt_data(:).sq_gest],[6 length(all_su_data)]);
after_lin_gests = reshape([after_filt_data(:).lin_gest],[6 length(all_su_data)]);
for ee = 1:n_expts
    ee
    su_set = find([all_su_data(:).expt_num] ==  str2num(expt_list{ee}(2:end)));
    sp = all_su_data(su_set(1)).after_mod.stim_params(1);
    nPix = sp.stim_dims(2);
    if ismember(expt_list{ee},{'G091','G092','G093'})
        dds = 0.67;
    else
        dds = 0.12;
    end
    rand_stim = zeros(n_frames,nPix);
    rand_mat = rand(n_frames,nPix);
    rand_stim(rand_mat <= dds/2) = -1;
    rand_stim(rand_mat >= 1-(dds/2)) = 1;
    stim_params = NMMcreate_stim_params([flen nPix],0.01);
    test_stim = create_time_embedding(rand_stim,stim_params);
    
    
    for ii = 1:length(su_set)
        Xtargs = [all_su_data(su_set(ii)).after_mod.mods(:).Xtarget];
        
        uncor_filters = [all_su_data(su_set(ii)).before_mod.mods(Xtargs==1).filtK];
        filt_out = test_stim*uncor_filters;
        filt_out(:,2:end) = filt_out(:,2:end).^2;
        filt_out_var = std(filt_out);
        filt_out_var = filt_out_var/nansum(filt_out_var);
        simp_before(su_set(ii)) = filt_out_var(1);
        
        cor_filters = [all_su_data(su_set(ii)).after_mod.mods(Xtargs==1).filtK];
        filt_out = test_stim*cor_filters;
        filt_out(:,2:end) = filt_out(:,2:end).^2;
        filt_out_var = std(filt_out);
        filt_out_var = filt_out_var/nansum(filt_out_var);
        simp(su_set(ii)) = filt_out_var(1);
        
        all_means = [after_filt_data(su_set(ii)).lin_mean after_filt_data(su_set(ii)).sq_mean ...
            after_filt_data(su_set(ii)).bsq_mean];
%         all_stds = [after_filt_data(su_set(ii)).lin_std after_filt_data(su_set(ii)).sq_std ...
%             after_filt_data(su_set(ii)).bsq_std];
        all_stds = [after_lin_gests(2,su_set(ii)) after_sq_gests(2,su_set(ii)) ...
            after_bsq_gests(2,su_set(ii))];
        cur_ford = after_filt_data(su_set(ii)).ford;
        weight_mean(su_set(ii)) = nansum(all_means.*filt_out_var(cur_ford));
        weight_std(su_set(ii)) = nansum(all_stds.*filt_out_var(cur_ford));
    end
end


%% PLOT BEFORE AFTER COMPARISONS
thresh_R2 = 0;

before_sq_gests = reshape([before_filt_data(used_sus).sq_gest],[6 length(used_sus)]);
before_lin_gests = reshape([before_filt_data(used_sus).lin_gest],[6 length(used_sus)]);
after_sq_gests = reshape([after_filt_data(used_sus).sq_gest],[6 length(used_sus)]);
after_lin_gests = reshape([after_filt_data(used_sus).lin_gest],[6 length(used_sus)]);

before_sq_gR2 = [before_filt_data(used_sus).sq_est_R2];
before_lin_gR2 = [before_filt_data(used_sus).lin_est_R2];
after_sq_gR2 = [after_filt_data(used_sus).sq_est_R2];
after_lin_gR2 = [after_filt_data(used_sus).lin_est_R2];

good_lin_filts = after_lin_gR2 >= thresh_R2;
good_sq_filts = after_sq_gR2 >= thresh_R2;
% good_lin_filts = after_lin_gR2 >= thresh_R2 & before_lin_gR2 >= thresh_R2;
% good_sq_filts = after_sq_gR2 >= thresh_R2 & before_sq_gR2 >= thresh_R2;

before_mod_R2 = [all_su_data(used_sus).before_R2];
after_mod_R2 = [all_su_data(used_sus).after_R2];

lin_amp_imp = after_lin_gests(5,good_lin_filts)./before_lin_gests(5,good_lin_filts);
sq_amp_imp = after_sq_gests(5,good_sq_filts)./before_sq_gests(5,good_sq_filts);

lin_width_imp = after_lin_gests(2,good_lin_filts)./before_lin_gests(2,good_lin_filts);
sq_width_imp = after_sq_gests(2,good_sq_filts)./before_sq_gests(2,good_sq_filts);

lin_SF_imp = after_lin_gests(3,good_lin_filts)./before_lin_gests(3,good_lin_filts);
sq_SF_imp = after_sq_gests(3,good_sq_filts)./before_sq_gests(3,good_sq_filts);


%%
close all
mSize = 12;

h1=figure();hold on
blin = robustfit(before_lin_gests(2,good_lin_filts),after_lin_gests(2,good_lin_filts),[],[],'off');
bsq = robustfit(before_sq_gests(2,good_sq_filts),after_sq_gests(2,good_sq_filts),[],[],'off');
plot(2*before_lin_gests(2,good_lin_filts),2*after_lin_gests(2,good_lin_filts),'.','markersize',mSize);
plot(2*before_sq_gests(2,good_sq_filts),2*after_sq_gests(2,good_sq_filts),'r.','markersize',mSize);
% legend({'Linear filters','Squared filters'},'Location','SouthEast');
% xlim([0 0.6]);
maxlim = max(max(xlim(),ylim()));
xlim([0 maxlim]); ylim([0 maxlim]);
line([0 maxlim],[0 maxlim],'color','k','linestyle','--');
xl = xlim(); xx = linspace(xl(1),xl(2),100);
plot(xx,blin*xx,'b--'); plot(xx,bsq*xx,'r--');
xlabel('RF width before (deg)'); ylabel('RF width after (deg)');


blin = robustfit(before_lin_gests(3,good_lin_filts),after_lin_gests(3,good_lin_filts),[],[],'off');
bsq = robustfit(before_sq_gests(3,good_sq_filts),after_sq_gests(3,good_sq_filts),[],[],'off');
h2=figure();hold on
plot(before_lin_gests(3,good_lin_filts),after_lin_gests(3,good_lin_filts),'.','markersize',mSize);
plot(before_sq_gests(3,good_sq_filts),after_sq_gests(3,good_sq_filts),'r.','markersize',mSize);
% legend({'Linear filters','Squared filters'},'Location','SouthEast');
maxlim = max(max(xlim(),ylim()));
xlim([0 maxlim]); ylim([0 maxlim]);
line([0 maxlim],[0 maxlim],'color','k','linestyle','--');
xl = xlim(); xx = linspace(xl(1),xl(2),100);
plot(xx,blin*xx,'b--'); plot(xx,bsq*xx,'r--');
xlabel('Spatial freq before (cyc/deg)'); ylabel('Spatial freq after (cyc/deg)');

blin = robustfit(before_lin_gests(5,good_lin_filts),after_lin_gests(5,good_lin_filts),[],[],'off');
bsq = robustfit(before_sq_gests(5,good_sq_filts),after_sq_gests(5,good_sq_filts),[],[],'off');
h3=figure();hold on
plot(before_lin_gests(5,good_lin_filts),after_lin_gests(5,good_lin_filts),'.','markersize',mSize);
plot(before_sq_gests(5,good_sq_filts),after_sq_gests(5,good_sq_filts),'r.','markersize',mSize);
% legend({'Linear filters','Squared filters'},'Location','SouthEast');
maxlim = max(max(xlim(),ylim()));
xlim([0 maxlim]); ylim([0 maxlim]);
line([0 maxlim],[0 maxlim],'color','k','linestyle','--');
xl = xlim(); xx = linspace(xl(1),xl(2),100);
plot(xx,blin*xx,'b--'); plot(xx,bsq*xx,'r--');
xlabel('Filter amplitude before'); ylabel('Filter amplitude after');

b = robustfit(before_mod_R2,after_mod_R2,[],[],'off');
h4 = figure(); hold on
plot(before_mod_R2,after_mod_R2,'k.','markersize',mSize);
maxlim = max(max(xlim(),ylim()));
xlim([0 maxlim]); ylim([0 maxlim]);
xl = xlim(); xx = linspace(xl(1),xl(2),100);
plot(xx,b*xx,'k--');
line([0 maxlim],[0 maxlim],'color','k','linestyle','--');
xlabel('R2 before'); ylabel('R2 after');

%% PRINT FILTER BEFORE/AFTER COMPARISONS
% fig_width = 3.27; 
fig_width = 2.75; 
rel_height = 0.85;

figufy(h1);
fname = [fig_dir 'RFsize_before_after2.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h2);
fname = [fig_dir 'pSF_before_after2.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h2);

figufy(h3);
fname = [fig_dir 'filtamp_before_after2.pdf'];
exportfig(h3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h3);

figufy(h4);
fname = [fig_dir 'modR2_before_after2.pdf'];
exportfig(h4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h4);

%%
b = robustfit(simp_before,simp,[],[],'off');
h5 = figure(); hold on
plot(simp_before,simp,'k.','markersize',mSize);
maxlim = max(max(xlim(),ylim()));
xlim([0 maxlim]); ylim([0 maxlim]);
xl = xlim(); xx = linspace(xl(1),xl(2),100);
plot(xx,b*xx,'k--');
line([0 maxlim],[0 maxlim],'color','k','linestyle','--');
xlabel('R2 before'); ylabel('R2 after');

%% RF width vs R2 improvement SQ FILTS
% rmpath('~/James_scripts/bruce/bruce_code/')

expt_nums = [all_su_data(used_sus).expt_num];

fig_width = 3.27; rel_height = 0.8;
mSize = 15;

lemExpts = [270 266 275 277];
% lemExptsf = [266 275 ];
% lemUnits = ismember(expt_nums,lemExpts) & after_sq_gR2 >= thresh_R2;
% lemUnitsf = ismember(expt_nums,lemExptsf) & after_sq_gR2 >= thresh_R2;
% jbeUnits = ~ismember(expt_nums,lemExpts) & after_sq_gR2 >= thresh_R2;
% allUnits = find(after_sq_gR2 >= thresh_R2);
lemUnits = ismember(expt_nums,lemExpts) ;
% lemUnitsf = ismember(expt_nums,lemExptsf);
jbeUnits = ~ismember(expt_nums,lemExpts) ;

% h = figure();
h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth2.fig');
hold on
% plot(after_sq_gests(2,jbeUnits)*2,R2_imp_fac(jbeUnits),'.','markersize',mSize*0.6);
% plot(after_sq_gests(2,lemUnits)*2,R2_imp_fac(lemUnits),'r.','markersize',mSize);
% plot(after_sq_gests(2,lemUnitsf)*2,R2_imp_fac(lemUnitsf),'g.','markersize',mSize);
% plot(weight_std(used_sus(jbeUnits))*2,R2_imp_fac(jbeUnits),'.','markersize',mSize*0.6);
% plot(weight_std(used_sus(lemUnits))*2,R2_imp_fac(lemUnits),'r.','markersize',mSize);
% plot(weight_std(used_sus(lemUnitsf))*2,R2_imp_fac(lemUnitsf),'g.','markersize',mSize);
% plot(weight_std(used_sus(jbeUnits))*2,R2_imp_fac_fix(jbeUnits),'.','markersize',mSize*0.6);
% plot(weight_std(used_sus(lemUnits))*2,R2_imp_fac_fix(lemUnits),'r.','markersize',mSize);
% plot(weight_std(used_sus(lemUnitsf))*2,R2_imp_fac_fix(lemUnitsf),'g.','markersize',mSize);
plot(weight_std(used_sus(jbeUnits))*2,R2_imp_fac_fix_noXV(jbeUnits),'.','markersize',mSize*0.6);
plot(weight_std(used_sus(lemUnits))*2,R2_imp_fac_fix_noXV(lemUnits),'r.','markersize',mSize);
figufy(h);
xlabel('RF width (deg)');
ylabel('R2 fold imp');
xlim([0.06 0.6]);
ylim([1 5]);

% F = @(x,xdata)x(1)*xdata.^(x(2));
% x0 = [4 -0.5];
% [x,resnorm] = lsqcurvefit(F,x0,after_sq_gests(2,allUnits)*2,R2_imp_fac(allUnits));
% xx = linspace(0.01,1,100);
% plot(xx,x(1)*xx.^x(2),'g','linewidth',2);
% xlim([0.04 0.6]);
set(gca,'xscale','log','yscale','log');

fname = [fig_dir 'R2imp_vs_width_SUfixonly_noXV.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% RF width vs eccentricity
load ~/Data/bruce/general_array_data/array_pos_data.mat
interp_ecc = sqrt(interp_x.^2 + interp_y.^2);

all_probe_nums = [all_su_data(used_sus).probe_num];
all_ecc = interp_ecc(all_probe_nums);
rf_eccs = [all_su_data(used_sus).rf_ecc];

lemExpts = [266 270 275 277];

allUnits = 1:length(used_sus);
lemUnits = ismember(expt_nums,lemExpts);
jbeUnits = ~ismember(expt_nums,lemExpts);
rf_eccs(jbeUnits) = all_ecc(jbeUnits);

jit_amp = 0.01;
X = rf_eccs + randn(size(rf_eccs))*jit_amp;
% Y = after_sq_gests(2,:)*2;
Y = weight_std*2;

h = figure(); hold on
plot(X(jbeUnits),Y(jbeUnits),'.');
plot(X(lemUnits),Y(lemUnits),'r.');

mdl = LinearModel.fit(X(allUnits),Y(allUnits),'linear','Robustopts','on');
xr = linspace(nanmin(X(allUnits)),nanmax(X(allUnits)),100);
[ypred,yci] = predict(mdl,xr(:));
plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
plot(xr,yci,'--','color',[0.2 0.8 0.2]);
xlabel('RF ecc (deg)');
ylabel('RF width (deg)');
set(gca,'yscale','log');
ylim([0.05 0.5]);

fig_width = 3.27; rel_height = 0.8;
figufy(h);
% fname = [fig_dir 'RFwidth_vs_ecc.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);
