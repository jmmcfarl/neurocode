clear all
close all

Expt_nums = [86];
for ee = 1:length(Expt_nums)
    Expt_num = Expt_nums(ee);
    Expt_name = sprintf('G%.3d',Expt_num);
    data_dir = ['~/Data/bruce/' Expt_name];
    cd(data_dir);
    
    anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
    hor_anal_name = 'monoc_eyecorr_hbar.mat';
    ver_anal_name = 'monoc_eyecorr_vbar.mat';
    hor_mods_name = 'monoc_eyecorr_hbar_mods.mat';
    ver_mods_name = 'monoc_eyecorr_vbar_mods.mat';
    
    cd(anal_dir)
    fprintf('\n Expt %d \n',Expt_num);
    % disp(exist(hor_anal_name,'file'))
    % disp(exist(ver_anal_name,'file'))
    
    %%
    load(hor_anal_name,'*R2*','*LLimp*','et_tr_set');
    load(hor_mods_name,'all_mod_SU*');
    
    [n_sus,n_drift_its,~] = size(dit_LLimp_LOO);
    tr_set = et_tr_set;
    su_inds = find(all_mod_SU(tr_set) > 0);
    for ii = 1:n_sus
        drift_LLimp_LOO(ii,:) = squeeze(dit_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
        fix_LLimp_L00(ii,:) = squeeze(it_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
        fix_LLimp_L00(ii,1) = squeeze(it_LLimp(1,tr_set(su_inds(ii))));
        drift_LLimp_LOO(ii,:) = squeeze(dit_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
        fix_LLimp_L00(ii,:) = squeeze(it_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
        fix_LLimp_L00(ii,1) = squeeze(it_LLimp(1,tr_set(su_inds(ii))));
    end
    fix_rel_LLimp_L00 = bsxfun(@rdivide,fix_LLimp_L00,fix_LLimp_L00(:,1));
    drift_rel_LLimp_L00 = bsxfun(@minus,drift_LLimp_LOO,drift_LLimp_LOO(:,1));
    drift_rel_LLimp_L00 = bsxfun(@rdivide,drift_rel_LLimp_L00,fix_LLimp_L00(:,1));
    drift_tot_LLimp_L00 = bsxfun(@rdivide,drift_LLimp_LOO,fix_LLimp_L00(:,1));
    
    hor_tot_imp = drift_tot_LLimp_L00(:,end);
    clear fix_* drift_*
   %%
   load(ver_anal_name,'*R2*','*LLimp*','et_tr_set');
   load(ver_mods_name,'all_mod_SU*');
   
   [n_sus,n_drift_its,~] = size(dit_LLimp_LOO);
   tr_set = et_tr_set;
   su_inds = find(all_mod_SU(tr_set) > 0);
   for ii = 1:n_sus
       drift_LLimp_LOO(ii,:) = squeeze(dit_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
       fix_LLimp_L00(ii,:) = squeeze(it_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
       fix_LLimp_L00(ii,1) = squeeze(it_LLimp(1,tr_set(su_inds(ii))));
       drift_LLimp_LOO(ii,:) = squeeze(dit_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
       fix_LLimp_L00(ii,:) = squeeze(it_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
       fix_LLimp_L00(ii,1) = squeeze(it_LLimp(1,tr_set(su_inds(ii))));
   end
   fix_rel_LLimp_L00 = bsxfun(@rdivide,fix_LLimp_L00,fix_LLimp_L00(:,1));
   drift_rel_LLimp_L00 = bsxfun(@minus,drift_LLimp_LOO,drift_LLimp_LOO(:,1));
   drift_rel_LLimp_L00 = bsxfun(@rdivide,drift_rel_LLimp_L00,fix_LLimp_L00(:,1));
   drift_tot_LLimp_L00 = bsxfun(@rdivide,drift_LLimp_LOO,fix_LLimp_L00(:,1));
   
   ver_tot_imp = drift_tot_LLimp_L00(:,end);
end

%%
Expt_nums = [275];
% for ee = 1:length(Expt_nums)
ee=1
    Expt_num = Expt_nums(ee);
    Expt_name = sprintf('M%.3d',Expt_num);
    data_dir = ['~/Data/bruce/' Expt_name];
    cd(data_dir);
    
    anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
    anal_name = 'monoc_eyecorr.mat';
    mods_name = 'monoc_eyecorr_mods.mat';
    
    cd(anal_dir)
    fprintf('\n Expt %d \n',Expt_num);
    % disp(exist(hor_anal_name,'file'))
    % disp(exist(ver_anal_name,'file'))
    
    %%
    load(anal_name,'*R2*','*LLimp*','et_tr_set');
    load(mods_name,'all_mod_SU*');
    
    [n_sus,n_drift_its,~] = size(dit_LLimp_LOO);
    tr_set = et_tr_set;
    su_inds = find(all_mod_SU(tr_set) > 0);
    for ii = 1:n_sus
        drift_LLimp_LOO(ii,:) = squeeze(dit_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
        fix_LLimp_L00(ii,:) = squeeze(it_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
        fix_LLimp_L00(ii,1) = squeeze(it_LLimp(1,tr_set(su_inds(ii))));
        drift_LLimp_LOO(ii,:) = squeeze(dit_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
        fix_LLimp_L00(ii,:) = squeeze(it_LLimp_LOO(ii,:,tr_set(su_inds(ii))));
        fix_LLimp_L00(ii,1) = squeeze(it_LLimp(1,tr_set(su_inds(ii))));
    end
    fix_rel_LLimp_L00 = bsxfun(@rdivide,fix_LLimp_L00,fix_LLimp_L00(:,1));
    drift_rel_LLimp_L00 = bsxfun(@minus,drift_LLimp_LOO,drift_LLimp_LOO(:,1));
    drift_rel_LLimp_L00 = bsxfun(@rdivide,drift_rel_LLimp_L00,fix_LLimp_L00(:,1));
    drift_tot_LLimp_L00 = bsxfun(@rdivide,drift_LLimp_LOO,fix_LLimp_L00(:,1));
    
    tot_imp = drift_tot_LLimp_L00(:,end);
    clear fix_* drift_*
% end