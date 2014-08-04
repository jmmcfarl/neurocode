   
clear all

Expt_name = 'G086';

data_dir = ['~/Data/bruce/' Expt_name];

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
cd(anal_dir);
fig_dir = '/home/james/Analysis/bruce/ET_final/';

bar_ori = 0;

if bar_ori == 0
    %for horizontal
    mod_data_name = 'monoc_eyecorr_hbar_mods.mat';
    anal_name = 'monoc_eyecorr_hbar.mat';
    mod_data_name_EC = 'monoc_eyecorr_hbar_mods_Cinit3.mat';
    mod_data_name_ECproc = 'monoc_eyecorr_hbar_mods_Cinit2.mat';
elseif bar_ori == 90
    %for vertical
    mod_data_name = 'monoc_eyecorr_vbar_mods.mat';
    anal_name = 'monoc_eyecorr_vbar.mat';
    mod_data_name_EC = 'monoc_eyecorr_vbar_mods_Cinit3.mat';
    mod_data_name_ECproc = 'monoc_eyecorr_vbar_mods_Cinit2.mat';
end

data = load(anal_name,'it_mods','it_R2','dit_mods_LOO','dit_mods','dit_R2_LOO','dit_R2','et_tr_set');
mod_data = load(mod_data_name,'all_mod_SU*');
mod_data_EC = load(mod_data_name_EC,'all_mod_R2','all_mod_fits');
mod_data_ECproc = load(mod_data_name_ECproc,'all_mod_R2','all_mod_fits');

clear cur_su_data
su_numbers = mod_data.all_mod_SUnum(data.et_tr_set(mod_data.all_mod_SUnum(data.et_tr_set) > 0));
su_inds = find(ismember(mod_data.all_mod_SUnum,su_numbers));
for ii = 1:length(su_numbers)
    cur_imp = data.dit_R2_LOO(ii,end,su_inds(ii));

    cur_su_data(ii).usable = true;
    cur_su_data(ii).before_mod = data.it_mods{1}(su_inds(ii));
    cur_su_data(ii).after_mod = data.dit_mods_LOO{ii,end}(su_inds(ii));
    cur_su_data(ii).before_R2 = data.it_R2(1,su_inds(ii));
    cur_su_data(ii).before_mod_EC = mod_data_EC.all_mod_fits(su_inds(ii));
    cur_su_data(ii).before_mod_ECproc = mod_data_ECproc.all_mod_fits(su_inds(ii));
    cur_su_data(ii).before_R2_EC = mod_data_EC.all_mod_R2(su_inds(ii));
    cur_su_data(ii).before_R2_ECproc = mod_data_ECproc.all_mod_R2(su_inds(ii));
    cur_su_data(ii).after_R2 = cur_imp;
    cur_su_data(ii).after_R2_noXV = data.dit_R2(end,su_inds(ii));
    cur_su_data(ii).probe_num = mod_data.all_mod_SU(su_inds(ii));
    cur_su_data(ii).dx = 0.0565/2;
    cur_su_data(ii).expt_num = Expt_name;
end

clear cur_all_data
for ii = 1:length(data.et_tr_set)
    cur_all_data(ii).before_mod = data.it_mods{1}(data.et_tr_set(ii));
    cur_all_data(ii).after_mod = data.dit_mods{end}(data.et_tr_set(ii));
    cur_all_data(ii).before_mod_EC = mod_data_EC.all_mod_fits(data.et_tr_set(ii));
    cur_all_data(ii).before_mod_ECproc = mod_data_ECproc.all_mod_fits(data.et_tr_set(ii));
    cur_all_data(ii).before_R2 = data.it_R2(1,data.et_tr_set(ii));
    cur_all_data(ii).before_R2_EC = mod_data_EC.all_mod_R2(data.et_tr_set(ii));
    cur_all_data(ii).before_R2_ECproc = mod_data_ECproc.all_mod_R2(data.et_tr_set(ii));
    cur_all_data(ii).after_R2_noXV = data.dit_R2(end,data.et_tr_set(ii));
    cur_all_data(ii).probe_num = mod_data.all_mod_SU(data.et_tr_set(ii));
    cur_all_data(ii).dx = 0.0565/2;
    cur_all_data(ii).expt_num = Expt_name;
end

%%
all_R2 = [[cur_su_data(:).before_R2]' [cur_su_data(:).before_R2_EC]' [cur_su_data(:).before_R2_ECproc]' [cur_su_data(:).after_R2]'];
all_R2 = bsxfun(@rdivide,all_R2(:,2:end),all_R2(:,1));

h = figure;
boxplot(all_R2,'labels',{'Coil corrected','Processed coil','Inferred'});
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
ylabel('Relative R2');

%%
fig_width = 3.27; rel_height = 0.8;
figufy(h);
fname = [fig_dir 'ver_modR2_init_compare_MU.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);

%%
all_R2 = [[cur_all_data(:).before_R2]' [cur_all_data(:).before_R2_EC]' [cur_all_data(:).before_R2_ECproc]' [cur_all_data(:).after_R2_noXV]'];
all_R2 = bsxfun(@rdivide,all_R2(:,2:end),all_R2(:,1));

h=figure;
boxplot(all_R2);
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
ylim([0 2.3]);
%%
unit_data = cur_su_data;
% unit_data = cur_all_data;
%%
close all
clear after_filt_data before_filt_data*
for ii = 1:length(unit_data)
    fprintf('Processing filters for unit %d of %d\n',ii,length(unit_data));
    
    Xtargs = [unit_data(ii).after_mod.mods(:).Xtarget];
    cor_filters = [unit_data(ii).after_mod.mods(Xtargs==1).filtK];
    cor_stim_params = unit_data(ii).after_mod.stim_params(1);
    after_filt_data(ii) = get_filter_properties(cor_filters,cor_stim_params,unit_data(ii).dx);
    
    Xtargs = [unit_data(ii).before_mod.mods(:).Xtarget];
    cor_filters = [unit_data(ii).before_mod.mods(Xtargs==1).filtK];
    cor_stim_params = unit_data(ii).before_mod.stim_params(1);
    before_filt_data(ii) = get_filter_properties(cor_filters,cor_stim_params,unit_data(ii).dx);
    
    Xtargs = [unit_data(ii).before_mod_EC.mods(:).Xtarget];
    cor_filters = [unit_data(ii).before_mod_EC.mods(Xtargs==1).filtK];
    cor_stim_params = unit_data(ii).before_mod_EC.stim_params(1);
    before_filt_data_EC(ii) = get_filter_properties(cor_filters,cor_stim_params,unit_data(ii).dx);
    
    Xtargs = [unit_data(ii).before_mod_ECproc.mods(:).Xtarget];
    cor_filters = [unit_data(ii).before_mod_ECproc.mods(Xtargs==1).filtK];
    cor_stim_params = unit_data(ii).before_mod_ECproc.stim_params(1);
    before_filt_data_ECproc(ii) = get_filter_properties(cor_filters,cor_stim_params,unit_data(ii).dx);

end


%%
before_sq_gests = reshape([before_filt_data(:).sq_gest],[6 length(unit_data)]);
before_lin_gests = reshape([before_filt_data(:).lin_gest],[6 length(unit_data)]);
beforeEC_sq_gests = reshape([before_filt_data_EC(:).sq_gest],[6 length(unit_data)]);
beforeEC_lin_gests = reshape([before_filt_data_EC(:).lin_gest],[6 length(unit_data)]);
beforeECproc_sq_gests = reshape([before_filt_data_ECproc(:).sq_gest],[6 length(unit_data)]);
beforeECproc_lin_gests = reshape([before_filt_data_ECproc(:).lin_gest],[6 length(unit_data)]);
after_sq_gests = reshape([after_filt_data(:).sq_gest],[6 length(unit_data)]);
after_lin_gests = reshape([after_filt_data(:).lin_gest],[6 length(unit_data)]);

after_sq_gR2 = [after_filt_data(:).sq_est_R2];
after_lin_gR2 = [after_filt_data(:).lin_est_R2];

thresh_R2 = 0.8;
good_lin_filts = after_lin_gR2 >= thresh_R2;
good_sq_filts = after_sq_gR2 >= thresh_R2;

all_lin_widths = [2*before_lin_gests(2,good_lin_filts)' 2*beforeEC_lin_gests(2,good_lin_filts)' 2*beforeECproc_lin_gests(2,good_lin_filts)' 2*after_lin_gests(2,good_lin_filts)'];
all_lin_widths = bsxfun(@rdivide,all_lin_widths(:,2:end),all_lin_widths(:,1));
all_sq_widths = [2*before_sq_gests(2,good_sq_filts)' 2*beforeEC_sq_gests(2,good_sq_filts)' 2*beforeECproc_sq_gests(2,good_sq_filts)' 2*after_sq_gests(2,good_sq_filts)'];
all_sq_widths = bsxfun(@rdivide,all_sq_widths(:,2:end),all_sq_widths(:,1));

all_lin_amps = [before_lin_gests(5,good_lin_filts)' beforeEC_lin_gests(5,good_lin_filts)' beforeECproc_lin_gests(5,good_lin_filts)' after_lin_gests(5,good_lin_filts)'];
all_lin_amps = bsxfun(@rdivide,all_lin_amps(:,2:end),all_lin_amps(:,1));
all_sq_amps = [before_sq_gests(5,good_sq_filts)' beforeEC_sq_gests(5,good_sq_filts)' beforeECproc_sq_gests(5,good_sq_filts)' after_sq_gests(5,good_sq_filts)'];
all_sq_amps = bsxfun(@rdivide,all_sq_amps(:,2:end),all_sq_amps(:,1));

%%
h = figure;
subplot(2,1,1);
boxplot(all_lin_widths,'labels',{'Coil corrected','Processed coil','Inferred'});
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
ylabel('Relative RF width');
subplot(2,1,2);
boxplot(all_sq_widths,'labels',{'Coil corrected','Processed coil','Inferred'});
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
ylabel('Relative RF width');

%%
fig_width = 3.27; rel_height = 1.6;
figufy(h);
fname = [fig_dir 'ver_RFwidth_init_compare.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);

%%

figure
subplot(2,1,1);
boxplot(all_lin_amps,'labels',{'Coil corrected','Processed coil','Inferred'});
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
ylabel('Relative filter amplitude');
subplot(2,1,2);
boxplot(all_sq_amps,'labels',{'Coil corrected','Processed coil','Inferred'});
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
ylabel('Relative filter amplitude');

%%
fig_width = 3.27; rel_height = 1.6;
figufy(h);
fname = [fig_dir 'ver_filtamp_init_compare.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);
