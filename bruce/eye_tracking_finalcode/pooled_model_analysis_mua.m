clear all
close all

addpath('~/James_scripts/bruce/eye_tracking_finalcode/');
% expt_list = {'G086'};
% expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
% expt_list = {'M232','M235','M239'};
% expt_list = {'M266','M270','M275','M277'};
expt_list = {'M266','M270','M275','M277','G085','G086','G087','G088','G089','G091','G093','G095'};
% expt_list = {'M232','M235','M239','M266','M270','M275','M277','G085','G086','G087','G088','G089','G091','G093','G095'};
print_out = false;

fig_dir = '/home/james/Analysis/bruce/ET_final/';

use_ori = [true true];
% use_ori = [true false];

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
       anal_name = 'monoc_eyecorr_hbar.mat';
       if exist(anal_name,'file') && use_ori(1)
           hor_data = load(anal_name,'it_mods','it_R2','dit_mods*','dit_R2*','et_tr_set');
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
       anal_name = 'monoc_eyecorr_vbar.mat';
       if exist(anal_name,'file') && use_ori(2)
           ver_data = load(anal_name,'it_mods','it_R2','dit_mods*','dit_R2*','et_tr_set');
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
       hor_mu_numbers = hor_data.et_tr_set;
       ver_mu_numbers = ver_data.et_tr_set;
       hor_su_numbers = hor_mod_data.all_mod_SUnum(hor_mu_numbers);
       ver_su_numbers = ver_mod_data.all_mod_SUnum(ver_mu_numbers);
       [tot_use_mus,IA] = unique([hor_mu_numbers ver_mu_numbers]);
       tot_use_sus = [hor_su_numbers ver_su_numbers];
       tot_use_sus = tot_use_sus(IA);
       for ii = 1:length(tot_use_mus)
           cur_hor_loc = hor_mu_numbers(hor_mu_numbers == tot_use_mus(ii));
           if isempty(cur_hor_loc);
               cur_hor_imp = -Inf;
           else
               cur_hor_imp = hor_data.dit_R2(end,cur_hor_loc);
           end;
           cur_ver_loc = ver_mu_numbers(ver_mu_numbers == tot_use_mus(ii));
           if isempty(cur_ver_loc)
               cur_ver_imp = -Inf;
           else
               cur_ver_imp = ver_data.dit_R2(end,cur_ver_loc);
           end
           if cur_hor_imp > cur_ver_imp
               cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 0;
               cur_su_data(ii).before_mod = hor_data.it_mods{1}(cur_hor_loc);
               cur_su_data(ii).after_mod = hor_data.dit_mods{end}(cur_hor_loc);
               cur_su_data(ii).before_R2 = hor_data.it_R2(1,cur_hor_loc);
               cur_su_data(ii).after_R2 = cur_hor_imp;
               cur_su_data(ii).after_R2_fix = hor_data.it_R2(end,cur_hor_loc);
               if cur_hor_loc > 96
                   cur_probe_num = hor_mod_data.all_mod_SU(cur_hor_loc);
                   cur_hor_spkloc = find(hor_SU_spk_numbers == tot_use_sus(ii));
                   cur_su_data(ii).num_blocks = cur_hor_SU_data.sua_data(cur_hor_spkloc).num_blocks;
                   cur_su_data(ii).n_spikes = cur_hor_SU_data.sua_data(cur_hor_spkloc).n_spikes;
                   cur_su_data(ii).mean_rate = cur_hor_SU_data.sua_data(cur_hor_spkloc).mean_rate;
               else
                   cur_probe_num = cur_hor_loc;
                   cur_su_data(ii).num_blocks = Inf;
                   cur_su_data(ii).n_spikes = Inf;
                   cur_su_data(ii).mean_rate = Inf;
               end
             cur_su_data(ii).probe_num = cur_probe_num;              
           elseif cur_ver_imp > cur_hor_imp
               cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 90;
               cur_su_data(ii).before_mod = ver_data.it_mods{1}(cur_ver_loc);
               cur_su_data(ii).after_mod = ver_data.dit_mods{end}(cur_ver_loc);
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
                   cur_su_data(ii).n_spikes = Inf;
                   cur_su_data(ii).mean_rate = Inf;
               end
               cur_su_data(ii).probe_num = cur_probe_num;                             
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
       anal_name = 'monoc_eyecorr_Cprior.mat';
       if exist(anal_name,'file')
           load(anal_name,'it_mods','it_R2','dit_mods','dit_R2','et_tr_set');
           load(mod_data_name,'all_mod_SU*');
       else
           error('No data');
       end
       clear cur_su_data
       mu_numbers = et_tr_set;
       su_numbers = all_mod_SUnum(mu_numbers);
       cur_SU_data = lem_SU_data{cellfun(@(x) x.Expt_num,lem_SU_data) == cur_expt_num};
       SU_spk_numbers = [cur_SU_data.sua_data(:).su_num];
      
       for ii = 1:length(mu_numbers)
           cur_loc = mu_numbers(ii);
           if ~isempty(cur_loc);
               cur_imp = dit_R2(end,cur_loc);
               cur_su_data(ii).usable = true;
               cur_su_data(ii).ori = 0;
               cur_su_data(ii).before_mod = it_mods{1}(cur_loc);
               cur_su_data(ii).after_mod = dit_mods{end}(cur_loc);
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
                  cur_su_data(ii).n_spikes = Inf;
                  cur_su_data(ii).mean_rate = Inf;
               end
               cur_su_data(ii).probe_num = cur_probe_num;
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
save pooled_mu_data all_su_data after_filt_data

%% COMPUTE SIMPLICITY
n_frames = 1e5;
flen = 12;
simp = nan(length(all_su_data),1);
weight_mean = nan(length(all_su_data),1);
weight_std = nan(length(all_su_data),1);
best_mean = nan(length(all_su_data),1);
best_std = nan(length(all_su_data),1);

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
        cor_filters = [all_su_data(su_set(ii)).after_mod.mods(Xtargs==1).filtK];
        filt_out = test_stim*cor_filters;
        filt_out(:,2:end) = filt_out(:,2:end).^2;
        filt_out_var = std(filt_out);
%         filt_out_var = var(filt_out);
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
        [~,best_filt] = max(filt_out_var);
        best_mean(su_set(ii)) = all_means(best_filt);
        best_std(su_set(ii)) = all_stds(best_filt);
    end
end

%%
min_num_blocks = 3;
num_blocks = [all_su_data(:).num_blocks];

used_units = find(num_blocks >= min_num_blocks);

after_R2 = [all_su_data(used_units).after_R2];
after_R2_fix = [all_su_data(used_units).after_R2_fix];
before_R2 = [all_su_data(used_units).before_R2];
R2_imp_fac = after_R2./before_R2;
R2_imp_fac_fix = after_R2_fix./before_R2;

%%
load ~/Data/bruce/general_array_data/array_pos_data.mat
interp_ecc = sqrt(interp_x.^2 + interp_y.^2);

all_probe_nums = [all_su_data(used_units).probe_num];
all_ecc = interp_ecc(all_probe_nums);

%%
expt_nums = [all_su_data(used_units).expt_num];
rf_eccs = [all_su_data(used_units).rf_ecc];
lemExpts = [266 270 275 277];
lemUnits = ismember(expt_nums,lemExpts);
jbeUnits = ~ismember(expt_nums,lemExpts);
rf_eccs(jbeUnits) = all_ecc(jbeUnits);

% [unique_expt_nums,IA,IC] = unique(expt_nums);
% for ii = 1:length(unique_expt_nums);
%     expt_rf_ecc(ii) = unique(rf_eccs(expt_nums==unique_expt_nums(ii)));
% end
% [~,expt_ord] = sort(expt_rf_ecc);
% G = expt_ord(IC);

after_sq_gests = reshape([after_filt_data(used_units).sq_gest],[6 length(used_units)]);
after_lin_gests = reshape([after_filt_data(used_units).lin_gest],[6 length(used_units)]);

thresh_R2 = 0.8;
after_sq_gR2 = [after_filt_data(used_units).sq_est_R2];
after_lin_gR2 = [after_filt_data(used_units).lin_est_R2];

%% RF width vs R2 improvement WEIGHTED WIDTH
% rmpath('~/James_scripts/bruce/bruce_code/')
 
fig_width = 3.27; rel_height = 0.8;
mSize = 10;

% allUnits = find(after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2);
allUnits = find(after_sq_gR2 >= thresh_R2 | after_lin_gR2 >= thresh_R2);
% allUnits = 1:length(used_units);

lemExpts = [266 270 275 277];
lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));

% h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth.fig');
h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth2.fig');
hold on
plot(weight_std(used_units(jbeUnits))*2,R2_imp_fac_fix(jbeUnits),'.','markersize',mSize*0.6);
plot(weight_std(used_units(lemUnits))*2,R2_imp_fac_fix(lemUnits),'r.','markersize',mSize);
line_hand = findall(gcf,'type','line','-and','color','k');
uistack(line_hand,'top');
figufy(h);

X = 2*weight_std(used_units);
% X = 2*after_sq_gests(2,:);
% Y = R2_imp_fac;
Y = R2_imp_fac_fix;
n_bins = 20;
pbin_edges = prctile(X,linspace(0,100,n_bins+1));
avg_X = nan(n_bins,1);
avg_Y = nan(n_bins,1);
std_X = nan(n_bins,1);
std_Y = nan(n_bins,1);
for ii = 1:n_bins
   curset = find(X >= pbin_edges(ii) & X < pbin_edges(ii+1));
   avg_X(ii) = mean(X(curset));
   avg_Y(ii) = mean(Y(curset));
   std_X(ii) = std(X(curset))/sqrt(length(curset));
   std_Y(ii) = std(Y(curset))/sqrt(length(curset));
%    std_X(ii) = std(X(curset));
%    std_Y(ii) = std(Y(curset));
end

% eh = errorbarxy(avg_X, avg_Y, std_X, std_Y,{'g' 'g' 'g'});
% set(eh(1),'linewidth',2);
% set(eh(2:end),'linewidth',2);

[~,xord] = sort(X(allUnits));
ysm = smooth(Y(allUnits(xord)),100,'moving');
plot(X(allUnits(xord)),ysm,'g-','linewidth',2);

xlim([0.06 0.6]);
ylim([1 5]);
set(gca,'xscale','log','yscale','log');

% F = @(x,xdata)x(1)*xdata.^(x(2));
% x0 = [4 -0.5];
% [x,resnorm] = lsqcurvefit(F,x0,X,Y');
% xx = linspace(0.01,1,100);
% plot(xx,x(1)*xx.^x(2),'g','linewidth',2);
% 
fname = [fig_dir 'modimp_vs_RFwidth_overlay3.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);



%% RF width vs R2 improvement SQ FILTS
% rmpath('~/James_scripts/bruce/bruce_code/')
 
fig_width = 3.27; rel_height = 0.8;
mSize = 10;

allUnits = find(after_sq_gR2 >= thresh_R2);
% lemExpts = [266 270 275 277];
% lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
% jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));
lemExptsf = [266 275];
lemExpts = [270 277];
lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
lemUnitsf = allUnits(ismember(expt_nums(allUnits),lemExptsf));
jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));

h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth.fig');
hold on
plot(after_sq_gests(2,jbeUnits)*2,R2_imp_fac(jbeUnits),'.','markersize',mSize*0.6);
plot(after_sq_gests(2,lemUnits)*2,R2_imp_fac(lemUnits),'r.','markersize',mSize);
plot(after_sq_gests(2,lemUnitsf)*2,R2_imp_fac(lemUnitsf),'g.','markersize',mSize);
line_hand = findall(gcf,'type','line','-and','color','k');
uistack(line_hand,'top');
figufy(h);

F = @(x,xdata)x(1)*xdata.^(x(2));
x0 = [4 -0.5];
% [x,resnorm] = lsqcurvefit(F,x0,after_sq_gests(2,allUnits)*2,R2_imp_fac(allUnits));
% xx = linspace(0.01,1,100);
% plot(xx,x(1)*xx.^x(2),'g','linewidth',2);
xlim([0.04 0.6]);
ylim([1 6.1])
set(gca,'xscale','log','yscale','log');

% fname = [fig_dir 'modimp_vs_RFwidth_overlay2.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);
% 
%% Lin RF width vs R2 improvement 
% rmpath('~/James_scripts/bruce/bruce_code/')
 
fig_width = 3.27; rel_height = 0.8;
mSize = 10;

lemExpts = [266 270 275 277];
lemUnits = ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
jbeUnits = ~ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
allUnits = find(after_lin_gR2 >= thresh_R2);

h=open('/home/james/Analysis/bruce/ET_final/modimp_vs_RFwidth.fig');
hold on
plot(after_lin_gests(2,jbeUnits)*2,R2_imp_fac(jbeUnits),'.','markersize',mSize*0.6);
plot(after_lin_gests(2,lemUnits)*2,R2_imp_fac(lemUnits),'r.','markersize',mSize);
line_hand = findall(gcf,'type','line','-and','color','k');
uistack(line_hand,'top');
figufy(h);

F = @(x,xdata)x(1)*xdata.^(x(2));
x0 = [8 -1.5];
[x,resnorm] = lsqcurvefit(F,x0,after_lin_gests(2,allUnits)*2,R2_imp_fac(allUnits));
xx = linspace(0.01,1,100);
plot(xx,x(1)*xx.^x(2),'r','linewidth',2);
xlim([0.04 0.8]);
set(gca,'xscale','log','yscale','log');

%% Lin SF vs R2 improvement 
fig_width = 3.27; rel_height = 0.8;
mSize = 10;

lemExpts = [266 270 275 277];
lemUnits = ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
jbeUnits = ~ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
allUnits = find(after_lin_gR2 >= thresh_R2);

h1 = figure(); hold on
plot(after_lin_gests(3,jbeUnits),R2_imp_fac(jbeUnits),'.','markersize',mSize*0.6);
plot(after_lin_gests(3,lemUnits),R2_imp_fac(lemUnits),'r.','markersize',mSize);
line_hand = findall(gcf,'type','line','-and','color','k');
uistack(line_hand,'top');

X = after_lin_gests(3,allUnits);
Y = R2_imp_fac(allUnits);
mdl = LinearModel.fit(X,Y,'linear','Robustopts','on');
xr = linspace(nanmin(X),nanmax(X),100);
[ypred,yci] = predict(mdl,xr(:));
plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
plot(xr,yci,'--','color',[0.2 0.8 0.2]);

figufy(h);
% fname = [fig_dir 'modimp_vs_LINRFwidth_overlay.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% Simplicity vs R2 improvement 
fig_width = 3.27; rel_height = 0.8;
mSize = 10;

lemExpts = [266 270 275 277];
lemUnits = ismember(expt_nums,lemExpts) & after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2;
jbeUnits = ~ismember(expt_nums,lemExpts) & after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2;
allUnits = after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2;

% X = log10(after_sq_gests(5,:)./after_lin_gests(5,:));
X = log10(after_sq_gests(5,:).^2./after_lin_gests(5,:));
Y = R2_imp_fac;

h1 = figure(); hold on
plot(X(jbeUnits),Y(jbeUnits),'.','markersize',mSize*0.6);
plot(X(lemUnits),Y(lemUnits),'r.','markersize',mSize);
% plot(atanh(after_sq_gests(5,jbeUnits)./after_lin_gests(5,jbeUnits)),R2_imp_fac(jbeUnits),'.','markersize',mSize*0.6);
% plot(atanh(after_sq_gests(5,lemUnits)./after_lin_gests(5,lemUnits)),R2_imp_fac(lemUnits),'r.','markersize',mSize);

% X = log10(after_sq_gests(5,allUnits)./after_lin_gests(5,allUnits));
X = log10(after_sq_gests(5,allUnits).^2./after_lin_gests(5,allUnits));
% X = after_lin_gests(5,allUnits);
Y = R2_imp_fac(allUnits);
mdl = LinearModel.fit(X,Y,'linear','Robustopts','on');
xr = linspace(nanmin(X),nanmax(X),100);
[ypred,yci] = predict(mdl,xr(:));
plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
plot(xr,yci,'--','color',[0.2 0.8 0.2]);
set(gca,'yscale','log')

fname = [fig_dir 'modimp_vs_simp.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);

%% Simplicity vs R2 improvement 
fig_width = 3.27; rel_height = 0.8;
mSize = 10;

lemExpts = [266 270 275 277];
lemUnits = ismember(expt_nums,lemExpts) & after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2;
jbeUnits = ~ismember(expt_nums,lemExpts) & after_sq_gR2 >= thresh_R2 & after_lin_gR2 >= thresh_R2;
% lemUnits = find(ismember(expt_nums,lemExpts));
% jbeUnits = find(~ismember(expt_nums,lemExpts));
allUnits = [jbeUnits lemUnits];

X = simp(used_units);
Y = R2_imp_fac;

h1 = figure(); hold on
plot(X(jbeUnits),Y(jbeUnits),'.','markersize',mSize*0.6);
plot(X(lemUnits),Y(lemUnits),'r.','markersize',mSize);

mdl = LinearModel.fit(X,Y,'linear','Robustopts','on');
xr = linspace(nanmin(X),nanmax(X),100);
[ypred,yci] = predict(mdl,xr(:));
plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
plot(xr,yci,'--','color',[0.2 0.8 0.2]);
yl = ylim();
ylim([1 6.1]);
set(gca,'yscale','log')

fname = [fig_dir 'modimp_vs_simp.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

%% RF width vs eccentricity
allUnits = find(after_sq_gR2 >= thresh_R2 | after_lin_gR2 >= thresh_R2);
% allUnits = 1:length(used_units);

lemExpts = [266 270 275 277];
% lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
% jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));

% allUnits = 1:length(used_units);
allUnits = find(after_sq_gR2 >= thresh_R2 | after_lin_gR2 >= thresh_R2);
% allUnits = find(after_sq_gR2 >= thresh_R2);
lemUnits = allUnits(ismember(expt_nums(allUnits),lemExpts));
jbeUnits = allUnits(~ismember(expt_nums(allUnits),lemExpts));

jit_amp = 0.01;
X = rf_eccs + randn(size(rf_eccs))*jit_amp;
% Y = after_sq_gests(2,:)*2;
Y = weight_std(used_units)*2;

h = figure(); hold on
plot(X(jbeUnits),Y(jbeUnits),'.');
plot(X(lemUnits),Y(lemUnits),'r.');

mdl = LinearModel.fit(X(allUnits),Y(allUnits),'linear','Robustopts','on');
% mdl = LinearModel.fit(X(allUnits),Y(allUnits),'linear','Robustopts','off');
xr = linspace(nanmin(X(allUnits)),nanmax(X(allUnits)),100);
[ypred,yci] = predict(mdl,xr(:));
plot(xr,ypred,'color',[0.2 0.8 0.2],'linewidth',2);
plot(xr,yci,'--','color',[0.2 0.8 0.2]);
xlabel('RF ecc (deg)');
ylabel('RF width (deg)');
% set(gca,'yscale','log');
ylim([0.05 0.5]);

fig_width = 3.27; rel_height = 0.8;
figufy(h);
fname = [fig_dir 'RFwidth_vs_ecc.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);

%% SF vs eccentricity
lemExpts = [266 270 275 277];

lemUnits = ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
jbeUnits = ~ismember(expt_nums,lemExpts) & after_lin_gR2 >= thresh_R2;
allUnits = find(after_lin_gR2 >= thresh_R2);

jit_amp = 0.01;
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
ylabel('RF width (deg)');
% set(gca,'yscale','log');
% ylim([0.05 0.5]);

fig_width = 3.27; rel_height = 0.8;
figufy(h);
% fname = [fig_dir 'RFwidth_vs_ecc.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

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

