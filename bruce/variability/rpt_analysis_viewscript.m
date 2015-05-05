clear all

% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
% expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
Expt_list = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};
expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
expt_mname = repmat({'lem'},1,length(Expt_list));
expt_rnum = ones(length(Expt_list),2);

Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012','M013'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 160; 0 0; 100 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe','jbe'});
expt_rnum = cat(1,expt_rnum,[1 1; 1 1; 1 1; 1 1; 1 2; 1 2; 1 1]);

%%
close all
base_sname = 'rpt_variability_compact';

Ccnt = 1;
Pcnt = 1;
all_Cdata = [];
all_Pdata = [];
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    monk_name = expt_mname{Elist_cnt};
    for bori_cnt = 1:2
        bar_ori = expt_oris(Elist_cnt,bori_cnt);
        rec_number = expt_rnum(Elist_cnt,bori_cnt);
        if ~isnan(bar_ori)
            fprintf('Loading %s on Expt %s ori %d\n',base_sname,Expt_name,bar_ori);
            data_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
            
            sname = [data_dir base_sname sprintf('_ori%d',bar_ori)];
            if rec_number > 1
               sname = strcat(sname,sprintf('_r%d',rec_number)); 
            end
            load(sname);
              
            if ~isfield(all_Cdata,'bestGQM') & isfield(EP_data,'bestGQM')
                EP_data = rmfield(EP_data,'bestGQM');
            end
            if ~isfield(all_Pdata,'EP_xcovar_LOO') & isfield(EP_pairs,'EP_xcovar_LOO')
                EP_pairs = rmfield(EP_pairs,'EP_xcovar_LOO');
            end
            
            for cc = 1:length(EP_data)
                if ~isnan(EP_data(cc).ov_avg_BS)
                    fprintf('Cell %d/%d\n',cc,length(EP_data));
                    
                    EP_data(cc).monkey = monk_name;
                    EP_data(cc).Expt_num = str2num(Expt_name(2:end));
                    EP_data(cc).bar_ori = bar_ori;
                     EP_data(cc).rec_number = rec_number;
                    EP_data(cc).cell_ID = Ccnt;
                    all_Cdata = cat(1,all_Cdata,EP_data(cc));
                    Ccnt = Ccnt + 1;                    
                end
            end
            
            cur_ucells = find(~isnan([EP_data(:).ov_avg_BS]));
            for cc = 1:length(EP_pairs)
                if all(ismember(EP_pairs(cc).ids,cur_ucells))
                    fprintf('Cell Pair %d/%d\n',cc,length(EP_pairs));
                    
                    EP_pairs(cc).monkey = monk_name;
                    EP_pairs(cc).Expt_num = str2num(Expt_name(2:end));
                    EP_pairs(cc).bar_ori = bar_ori;
                    EP_pairs(cc).rec_number = rec_number;
                    EP_pairs(cc).cell_IDs = [EP_data([EP_pairs(cc).ids]).cell_ID];
                    all_Pdata = cat(1,all_Pdata,EP_pairs(cc));
                    Pcnt = Pcnt + 1;
                end
            end
        end
    end
end

%% FOR CELLS RECORDED MULTIPLE TIMES, PICK BEST INSTANCE
SU_numbers = arrayfun(@(x) x.unit_data.SU_number,all_Cdata);
Expt_numbers = [all_Cdata(:).Expt_num]';
Rec_numbers = [all_Cdata(:).rec_number]';
to_eliminate = [];
for ii = 1:length(all_Cdata)
    curset = find(SU_numbers == SU_numbers(ii) & Expt_numbers == Expt_numbers(ii) & Rec_numbers == Rec_numbers(ii));
    if length(curset) > 1
%         cur_xvLLs = arrayfun(@(x) x.ModData.bestGQM.xvLLimp,all_Cdata(curset));
        avg_rates = arrayfun(@(x) x.unit_data.avg_rate,all_Cdata(curset));
%         xvLL_rate = cur_xvLLs.*avg_rates;
        
        [~,best_ind] = max(avg_rates);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
    end
end

to_eliminate = unique(to_eliminate);
double_CIDs = [all_Cdata(to_eliminate).cell_ID];
fprintf('Eliminating %d/%d duplicate SUs (multiple oris)\n',length(to_eliminate),length(all_Cdata));
all_Cdata(to_eliminate) = [];
pair_IDs = cat(1,all_Pdata.cell_IDs);
all_Pdata(any(ismember(pair_IDs,to_eliminate),2)) = [];

%% FOR SAME SUS RECORDED ON MULTIPLE SESSIONS WITH DIFFERENT ED
dup_SUs = [12 1 5; 12 3 8]; %[Expt_num r2_SU_Number r1_SU_number]
to_eliminate = [];
for ii = 1:size(dup_SUs,1)
   cur_unit_1 = find(Expt_numbers == dup_SUs(ii,1) & Rec_numbers == 2 & SU_numbers == dup_SUs(ii,2));
   cur_unit_2 = find(Expt_numbers == dup_SUs(ii,1) & Rec_numbers == 1 & SU_numbers == dup_SUs(ii,3));
   curset = [cur_unit_1 cur_unit_2];
   if length(curset) == 2
%         cur_xvLLs = arrayfun(@(x) x.ModData.bestGQM.xvLLimp,all_Cdata(curset));
        avg_rates = arrayfun(@(x) x.unit_data.avg_rate,all_Cdata(curset));
%         xvLL_rate = cur_xvLLs.*avg_rates;
        [~,best_ind] = max(avg_rates);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
   end
end

double_CIDs = [all_Cdata(to_eliminate).cell_ID];
fprintf('Eliminating %d/%d duplicate SUs (multiple recs)\n',length(to_eliminate),length(all_Cdata));
elim_CIDs = [all_Cdata(to_eliminate).cell_ID];
all_Cdata(to_eliminate) = [];
pair_IDs = cat(1,all_Pdata.cell_IDs);
all_Pdata(any(ismember(pair_IDs,elim_CIDs),2)) = [];

%%

all_ntrials = arrayfun(@(x) sum(x.n_utrials),all_Cdata);
all_avgrates = [all_Cdata(:).ov_avg_BS]'/EP_params.base_dt;
all_spline_vars = arrayfun(@(x) x.spline_pred_looEP(1),all_Cdata);
all_spline_vars_noLOO = arrayfun(@(x) x.spline_pred_baseEP(1),all_Cdata);

all_psth_vars = [all_Cdata(:).pair_psth_var]';
% all_ball_vars = cat(2,all_Cdata.eps_ball_var)';
% all_ball_alphas = bsxfun(@rdivide,all_psth_vars,all_ball_vars);
all_spline_alphas = all_psth_vars./all_spline_vars;
all_spline_alphas_noLOO = all_psth_vars./all_spline_vars_noLOO;

all_mod_psth_vars = arrayfun(@(x) mean(x.mod_psth_vars),all_Cdata);
all_mod_tot_vars = arrayfun(@(x) mean(x.mod_tot_vars),all_Cdata);
all_mod_alphas = all_mod_psth_vars./all_mod_tot_vars;

all_psth_FF = [all_Cdata(:).psth_FF];
all_spline_FF = [all_Cdata(:).spline_FF];

all_monkey = {all_Cdata(:).monkey};
all_CID = [all_Cdata(:).cell_ID];

pair_IDs = cat(1,all_Pdata.cell_IDs);
pair_matches = nan(size(pair_IDs));
for ii = 1:size(pair_IDs,1)
    pair_matches(ii,1) = find(all_CID==pair_IDs(ii,1));
    pair_matches(ii,2) = find(all_CID == pair_IDs(ii,2));
end


RF_ecc = arrayfun(@(x) x.tune_props.RF_ecc,all_Cdata);
RF_width = 2*arrayfun(@(x) x.tune_props.RF_sigma,all_Cdata);
RF_PSF = arrayfun(@(x) x.tune_props.RF_FSF,all_Cdata);
RF_PRM = arrayfun(@(x) x.tune_props.PRM,all_Cdata);

min_nTrials = 50;
min_avgRate = 5;
uset = find(all_ntrials >= min_nTrials & all_avgrates >= min_avgRate);
upairs = find(all(ismember(pair_IDs,all_CID(uset)),2) & pair_IDs(:,1) ~= pair_IDs(:,2));

%%
% close all

f1 = figure(); hold on
jbe_units = uset(strcmp(all_monkey(uset),'jbe'));
lem_units = uset(strcmp(all_monkey(uset),'lem'));
plot(all_mod_alphas(lem_units),all_spline_alphas(lem_units),'o');
plot(all_mod_alphas(jbe_units),all_spline_alphas(jbe_units),'ro');

% plot(all_mod_alphas(lem_units),all_ball_alphas(lem_units,4),'o');
% plot(all_mod_alphas(jbe_units),all_ball_alphas(jbe_units,4),'ro');
line([0 1],[0 1],'color','k');

%%
f2 = figure();
subplot(2,2,1)
plot(RF_ecc(uset),all_spline_alphas(uset),'o')
subplot(2,2,2)
plot(RF_width(uset),all_spline_alphas(uset),'o')
set(gca,'xscale','log'); xlim([0.05 2])
subplot(2,2,3)
plot(RF_PSF(uset),all_spline_alphas(uset),'o')
subplot(2,2,4)
plot(RF_PRM(uset),all_spline_alphas(uset),'o')

%%
all_Cvars = arrayfun(@(x) mean(x.tot_var),all_Cdata);

all_tot_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.tot_xcovar,1),all_Pdata(upairs),'uniformoutput',0)));
bad_pairs = find(isnan(all_tot_xcovs(:,11)));
upairs(bad_pairs) = []; all_tot_xcovs(bad_pairs,:) = [];

all_psth_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.pair_xcovar,1),all_Pdata(upairs),'uniformoutput',0)));
all_EP_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.EP_xcovar,1),all_Pdata(upairs),'uniformoutput',0)));

all_norms = sqrt(all_Cvars(pair_matches(upairs,1)).*all_Cvars(pair_matches(upairs,2)));

all_psth_noisecorr = bsxfun(@rdivide,all_tot_xcovs - all_psth_xcovs,all_norms);
all_EP_noisecorr = bsxfun(@rdivide,all_tot_xcovs - all_EP_xcovs,all_norms);
all_psth_sigcorr = bsxfun(@rdivide,all_psth_xcovs,all_norms);
all_EP_sigcorr = bsxfun(@rdivide,all_EP_xcovs,all_norms);

f1 = figure(); 
subplot(2,2,1);hold on
plot(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11),'.')
r1 = robustfit(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11));
xx = linspace(-0.2,0.2,100);
plot(xx,r1(1) + r1(2)*xx,'k')
 xlim([-0.3 0.3]); ylim([-0.3 0.3]);
 line([-0.3 0.3],[0 0],'color','k'); line([0 0],[-0.3 0.3],'color','k');
subplot(2,2,2);hold on
 plot(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11),'r.')
r2 = robustfit(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11));
 plot(xx,r2(1) + r2(2)*xx,'g')
xlim([-0.3 0.3]); ylim([-0.3 0.3]);
 line([-0.3 0.3],[0 0],'color','k'); line([0 0],[-0.3 0.3],'color','k');
subplot(2,2,3);hold on
plot(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11),'.')
r1 = robustfit(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11));
xx = linspace(-0.2,0.2,100);
plot(xx,r1(1) + r1(2)*xx,'k')
xlim([-0.15 0.15]); ylim([-0.15 0.15]);
 line([-0.3 0.3],[0 0],'color','k'); line([0 0],[-0.3 0.3],'color','k');
subplot(2,2,4);hold on
 plot(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11),'r.')
r2 = robustfit(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11));
 plot(xx,r2(1) + r2(2)*xx,'g')
xlim([-0.15 0.15]); ylim([-0.15 1.15]);
 line([-0.3 0.3],[0 0],'color','k'); line([0 0],[-0.3 0.3],'color','k');
%%

% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
% expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
Expt_list = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};
expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
expt_mname = repmat({'lem'},1,length(Expt_list));

Expt_list = cat(2,Expt_list,{'G085','G086','G087','G088','G089','G091','G093','G095'});
expt_oris = cat(1,expt_oris,[0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan]);
expt_mname = cat(2,expt_mname,repmat({'jbe'},1,8));

Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 nan; 0 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe'});

base_sname = 'model_variability_compact';

Mcnt = 1;
all_Mdata = [];
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    monk_name = expt_mname{Elist_cnt};
    for bori_cnt = 1:2
        bar_ori = expt_oris(Elist_cnt,bori_cnt);
        if ~isnan(bar_ori)
            fprintf('Loading %s on Expt %s ori %d\n',base_sname,Expt_name,bar_ori);
            data_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
            
            sname = [data_dir base_sname sprintf('_ori%d.mat',bar_ori)];
            load(sname);
                        
            EP_data = EP_data(targs);
            for cc = 1:length(EP_data)
                if ~isempty(EP_data(cc).ModData.unit_data)
                    fprintf('Cell %d/%d\n',cc,length(EP_data));
                    
                    EP_data(cc).monkey = monk_name;
                    EP_data(cc).Expt_num = str2num(Expt_name(2:end));
                    EP_data(cc).bar_ori = bar_ori;
                    EP_data(cc).cell_ID = Ccnt;
                    EP_data(cc).ov_EP_xcov = ov_EP_data.EP_xcov;
                    EP_data(cc).ov_EP_lags = ov_EP_data.EP_lags;
                    all_Mdata = cat(1,all_Mdata,EP_data(cc));
                    Mcnt = Mcnt + 1;                    
                end
            end
            
        end
    end
end

%%
SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,all_Mdata);
Expt_numbers = [all_Mdata(:).Expt_num]';
to_eliminate = [];
for ii = 1:length(all_Mdata)
    curset = find(SU_numbers == SU_numbers(ii) & Expt_numbers == Expt_numbers(ii));
    if length(curset) > 1
        cur_xvLLs = arrayfun(@(x) x.ModData.bestGQM.xvLLimp,all_Mdata(curset));
        avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_Mdata(curset));
        xvLL_rate = cur_xvLLs.*avg_rates;
        
        [~,best_ind] = max(cur_xvLLs);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
    end
end
to_eliminate = unique(to_eliminate);
fprintf('Eliminating %d/%d duplicate SUs\n',length(to_eliminate),length(all_Mdata));
all_Mdata(to_eliminate) = [];

%%
all_avgrates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_Mdata);
all_monkey = {all_Mdata(:).monkey};
all_CID = [all_Mdata(:).cell_ID];

RF_ecc = arrayfun(@(x) x.ModData.tune_props.RF_ecc,all_Mdata);
RF_width = 2*arrayfun(@(x) x.ModData.tune_props.RF_sigma,all_Mdata);
RF_PSF = arrayfun(@(x) x.ModData.tune_props.RF_FSF,all_Mdata);
RF_PRM = arrayfun(@(x) x.ModData.tune_props.PRM,all_Mdata);

actual_EP_SDs = arrayfun(@(x) x.poss_SDs(end),all_Mdata);

mod_alpha_funs = cat(1,all_Mdata.alpha_funs);

min_avgRate = 5;
% uset = find(all_avgrates >= min_avgRate & RF_ecc > 1);
uset = find(all_avgrates >= min_avgRate);

%%
ov_EP_xcov = cat(2,all_Mdata(uset).ov_EP_xcov);
EP_lags = all_Mdata(1).ov_EP_lags;
figure;
shadedErrorBar(EP_lags*.01,nanmean(ov_EP_xcov,2),nanstd(ov_EP_xcov,[],2));
%%
ep_FF_funs = cat(3,all_Mdata.ep_FF_est);
poss_ubins = all_Mdata(1).poss_ubins;
rate_covs = arrayfun(@(x) x.base_vars(end),all_Mdata);
rate_FF = rate_covs./all_avgrates;

f1 = figure();
errorbar(poss_ubins*0.01,squeeze(nanmean(ep_FF_funs(:,end,uset),3)),squeeze(nanstd(ep_FF_funs(:,end,uset),[],3))/sqrt(length(uset)))

%%
f2 = figure();
poss_SDs = all_Mdata(1).poss_SDs(1:end-1);
errorbar(poss_SDs,nanmean(mod_alpha_funs(uset,1:end-1)),nanstd(mod_alpha_funs(uset,1:end-1))/sqrt(length(uset)));
