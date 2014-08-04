clear all
close all

addpath('~/James_scripts/bruce/eye_tracking_finalcode/');
expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};

fig_dir = '/home/james/Analysis/bruce/ET_final/';

%%
n_expts = length(expt_list);
clear tot_LL_*
for ee = 1:n_expts
    ee
    cur_expt_name = expt_list{ee};
    cur_expt_num = str2num(cur_expt_name(2:end));
        
    data_dir = ['~/Data/bruce/' cur_expt_name];
        load([data_dir sprintf('/jbe%sExpts.mat',cur_expt_name)]); %load in Expts struct
       
       anal_dir = ['~/Analysis/bruce/' cur_expt_name '/ET_final/'];
       cd(anal_dir);
       
        clear dit_*
      %for horizontal
       mod_data_name = 'monoc_eyecorr_hbar_mods.mat';
%        anal_name = 'monoc_eyecorr_hbar.mat';
       anal_name = 'monoc_eyecorr_hbar.mat';
           load(anal_name,'it_mods','it_LL*','dit_LL*','et_tr_set');
       
        load('monoc_eyecorr_hbar_nspks');
        
        it_LL_tot = bsxfun(@times,it_LLimp(:,et_tr_set),tot_nspks);
        dit_LL_tot = bsxfun(@times,dit_LLimp(:,et_tr_set),tot_nspks);
        
        tot_LL_imp = [nansum(it_LL_tot,2); nansum(dit_LL_tot(2:end,:),2)];
        tot_LL_imp = tot_LL_imp - tot_LL_imp(1);
        tot_LL_imp_perspk(ee,:) = tot_LL_imp/sum(tot_nspks);
        
        nspk_frac = tot_nspks/sum(tot_nspks);
        tot_LL_imp_weight(ee,:) = [sum(bsxfun(@times,it_LLimp(:,et_tr_set),nspk_frac),2); sum(bsxfun(@times,dit_LLimp(:,et_tr_set),nspk_frac),2)];
        tot_LL_imp_weight(ee,:) = tot_LL_imp_weight(ee,:) - tot_LL_imp_weight(ee,1);
end

%%
fig_width = 4; rel_height = 0.8;
h =figure;
errorbar((1:size(tot_LL_imp_perspk,2))-1,mean(tot_LL_imp_perspk),std(tot_LL_imp_perspk),'ko-','markersize',8);
xlabel('Iteration number');
ylabel('LL improvement bits/spk')
figufy(h);
fname = [fig_dir 'EM_convergence.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);
