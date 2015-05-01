clear all

% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
% expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
Expt_list = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};
expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
expt_mname = repmat({'lem'},1,length(Expt_list));

Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 nan; 0 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe'});

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
        if ~isnan(bar_ori)
            fprintf('Loading %s on Expt %s ori %d\n',base_sname,Expt_name,bar_ori);
            data_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
            
            sname = [data_dir base_sname sprintf('_ori%d.mat',bar_ori)];
            load(sname);
                        
            for cc = 1:length(EP_data)
                if ~isnan(EP_data(cc).ov_avg_BS)
                    fprintf('Cell %d/%d\n',cc,length(EP_data));
                    
                    EP_data(cc).monkey = monk_name;
                    EP_data(cc).Expt_num = str2num(Expt_name(2:end));
                    EP_data(cc).bar_ori = bar_ori;
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
                    EP_pairs(cc).cell_IDs = [EP_data([EP_pairs(cc).ids]).cell_ID];
                    all_Pdata = cat(1,all_Pdata,EP_pairs(cc));
                    Pcnt = Pcnt + 1;
                end
            end
        end
    end
end

%% FOR CELLS RECORDED MULTIPLE TIMES, PICK BEST INSTANCE
SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,all_Cdata);
Expt_numbers = [all_Cdata(:).Expt_num]';
to_eliminate = [];
for ii = 1:length(all_Cdata)
    curset = find(SU_numbers == SU_numbers(ii) & Expt_numbers == Expt_numbers(ii));
    if length(curset) > 1
        cur_xvLLs = arrayfun(@(x) x.ModData.bestGQM.xvLLimp,all_Cdata(curset));
        avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_Cdata(curset));
        xvLL_rate = cur_xvLLs.*avg_rates;
        
        [~,best_ind] = max(cur_xvLLs);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
    end
end

to_eliminate = unique(to_eliminate);
double_CIDs = [all_Cdata(to_eliminate).cell_ID];
fprintf('Eliminating %d/%d duplicate SUs\n',length(to_eliminate),length(all_Mdata));
all_Cdata(to_eliminate) = [];

%%

all_ntrials = arrayfun(@(x) sum(x.n_utrials),all_Cdata);
all_avgrates = [all_Cdata(:).ov_avg_BS]'/EP_params.base_dt;
all_spline_vars = arrayfun(@(x) x.spline_pred_looEP(1),all_Cdata);

all_psth_vars = [all_Cdata(:).pair_psth_var]';
all_ball_vars = cat(2,all_Cdata.eps_ball_var)';
all_ball_alphas = bsxfun(@rdivide,all_psth_vars,all_ball_vars);
all_spline_alphas = all_psth_vars./all_spline_vars;

all_mod_psth_vars = arrayfun(@(x) mean(x.mod_psth_vars),all_Cdata);
all_mod_tot_vars = arrayfun(@(x) mean(x.mod_tot_vars),all_Cdata);
all_mod_alphas = all_mod_psth_vars./all_mod_tot_vars;

all_psth_FF = [all_Cdata(:).psth_FF];
all_spline_FF = [all_Cdata(:).spline_FF];

all_monkey = {all_Cdata(:).monkey};
all_CID = [all_Cdata(:).cell_ID];

pair_IDs = cat(1,all_Pdata.cell_IDs);



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
set(gca,'xscale','log'); xlim([0.05 1])
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


all_norms = sqrt(all_Cvars(pair_IDs(upairs,1)).*all_Cvars(pair_IDs(upairs,2)));

all_psth_noisecorr = bsxfun(@rdivide,all_tot_xcovs - all_psth_xcovs,all_norms);
all_EP_noisecorr = bsxfun(@rdivide,all_tot_xcovs - all_EP_xcovs,all_norms);
all_psth_sigcorr = bsxfun(@rdivide,all_psth_xcovs,all_norms);
all_EP_sigcorr = bsxfun(@rdivide,all_EP_xcovs,all_norms);

f1 = figure(); hold on
plot(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11),'.')
plot(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11),'r.')
r1 = robustfit(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11));
r2 = robustfit(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11));
xx = linspace(-0.2,0.2,100);
plot(xx,r1(1) + r1(2)*xx,'k')
 plot(xx,r2(1) + r2(2)*xx,'g')
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

mod_alpha_funs = cat(1,all_Mdata.alpha_funs);

min_avgRate = 5;
uset = find(all_avgrates >= min_avgRate);

