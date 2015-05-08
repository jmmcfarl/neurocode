clear all

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
% Expt_list = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};
% expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
expt_mname = repmat({'lem'},1,length(Expt_list));
expt_rnum = ones(length(Expt_list),2);

Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012','M013','M014'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 160; 0 0; 100 nan; 40 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe','jbe','jbe'});
expt_rnum = cat(1,expt_rnum,[1 1; 1 1; 1 1; 1 1; 1 2; 1 2; 1 1; 1 1]);

fig_dir = '/home/james/Analysis/bruce/variability/figures/';

%% load repeat trial data
close all
% base_sname = 'rpt_variability_compact';
base_sname = 'rpt_variability_compact_multDT';

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
                                 
            for cc = 1:size(EP_data,1)
                if ~isnan(EP_data(cc,1).ov_avg_BS)
                    fprintf('Cell %d/%d\n',cc,size(EP_data,1));
                    
                    EP_data(cc,1).monkey = monk_name;
                    EP_data(cc,1).Expt_num = str2num(Expt_name(2:end));
                    EP_data(cc,1).bar_ori = bar_ori;
                     EP_data(cc,1).rec_number = rec_number;
                    EP_data(cc,1).cell_ID = Ccnt;
                    all_Cdata = cat(1,all_Cdata,EP_data(cc,:));
                    Ccnt = Ccnt + 1;                    
                end
            end
            
            cur_ucells = find(~isnan([EP_data(:,1).ov_avg_BS]));
            for cc = 1:size(EP_pairs,1)
                if all(ismember(EP_pairs(cc,1).ids,cur_ucells))
%                     fprintf('Cell Pair %d/%d\n',cc,length(EP_pairs));
                    
                    EP_pairs(cc,1).monkey = monk_name;
                    EP_pairs(cc,1).Expt_num = str2num(Expt_name(2:end));
                    EP_pairs(cc,1).bar_ori = bar_ori;
                    EP_pairs(cc,1).rec_number = rec_number;
                    EP_pairs(cc,1).cell_IDs = [EP_data([EP_pairs(cc,1).ids]).cell_ID];
                    all_Pdata = cat(1,all_Pdata,EP_pairs(cc,:));
                    Pcnt = Pcnt + 1;
                end
            end
        end
    end
end
%%

% FOR CELLS RECORDED MULTIPLE TIMES, PICK BEST INSTANCE
SU_numbers = arrayfun(@(x) x.unit_data.SU_number,all_Cdata(:,1));

Expt_numbers = [all_Cdata(:,1).Expt_num]';
Rec_numbers = [all_Cdata(:,1).rec_number]';
to_eliminate = [];
for ii = 1:size(all_Cdata,1)
    curset = find(SU_numbers == SU_numbers(ii) & Expt_numbers == Expt_numbers(ii) & Rec_numbers == Rec_numbers(ii));
    if length(curset) > 1
        cur_xvLLs = arrayfun(@(x) x.bestGQM.xvLLimp,all_Cdata(curset,1));
        null_xvLLs = arrayfun(@(x) x.nullMod.xvLLimp,all_Cdata(curset,1));
        cur_xvLLs = cur_xvLLs - null_xvLLs;
        
        avg_rates = arrayfun(@(x) x.unit_data.avg_rate,all_Cdata(curset,1));
        xvLL_rate = cur_xvLLs.*avg_rates;
        
%         [~,best_ind] = max(avg_rates);
        [~,best_ind] = max(xvLL_rate);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
    end
end

to_eliminate = unique(to_eliminate);
fprintf('Eliminating %d/%d duplicate SUs (multiple oris)\n',length(to_eliminate),size(all_Cdata,1));
all_Cdata(to_eliminate,:) = [];
pair_IDs = cat(1,all_Pdata(:,1).cell_IDs);
all_Pdata(any(ismember(pair_IDs,to_eliminate),2),:) = [];

% FOR SAME SUS RECORDED ON MULTIPLE SESSIONS WITH DIFFERENT ED
dup_SUs = [12 1 5; 12 3 8]; %[Expt_num r2_SU_Number r1_SU_number]
SU_numbers = arrayfun(@(x) x.unit_data.SU_number,all_Cdata(:,1));
Expt_numbers = [all_Cdata(:,1).Expt_num]';
Rec_numbers = [all_Cdata(:,1).rec_number]';

to_eliminate = [];
for ii = 1:size(dup_SUs,1)
   cur_unit_1 = find(Expt_numbers == dup_SUs(ii,1) & Rec_numbers == 2 & SU_numbers == dup_SUs(ii,2));
   cur_unit_2 = find(Expt_numbers == dup_SUs(ii,1) & Rec_numbers == 1 & SU_numbers == dup_SUs(ii,3));
   curset = [cur_unit_1 cur_unit_2];
   if length(curset) == 2
        cur_xvLLs = arrayfun(@(x) x.bestGQM.xvLLimp,all_Cdata(curset,1));
        null_xvLLs = arrayfun(@(x) x.nullMod.xvLLimp,all_Cdata(curset,1));
        cur_xvLLs = cur_xvLLs - null_xvLLs;
        avg_rates = arrayfun(@(x) x.unit_data.avg_rate,all_Cdata(curset,1));
        xvLL_rate = cur_xvLLs.*avg_rates;
%         [~,best_ind] = max(avg_rates);
        [~,best_ind] = max(xvLL_rate);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
   end
end

fprintf('Eliminating %d/%d duplicate SUs (multiple recs)\n',length(to_eliminate),size(all_Cdata,1));
elim_CIDs = [all_Cdata(to_eliminate,1).cell_ID];
all_Cdata(to_eliminate,:) = [];
pair_IDs = cat(1,all_Pdata(:,1).cell_IDs);
all_Pdata(any(ismember(pair_IDs,elim_CIDs),2),:) = [];

%% extract SU properties and select units for analysis
poss_bin_dts = EP_params.poss_bin_dts;
poss_eps_sizes = EP_params.poss_eps_sizes;
n_SUs = size(all_Cdata,1);

all_ntrials = arrayfun(@(x) sum(x.n_utrials),all_Cdata(:,1));
all_avgrates = [all_Cdata(:,1).ov_avg_BS]'/poss_bin_dts(1);
all_spline_vars = arrayfun(@(x) x.spline_pred_looEP(1),all_Cdata);
all_spline_vars_noLOO = arrayfun(@(x) x.spline_pred_baseEP(1),all_Cdata);

all_psth_vars = arrayfun(@(x) x.pair_psth_var,all_Cdata);
all_bar_vars = nan(n_SUs,length(poss_bin_dts),length(poss_eps_sizes));
for ee = 1:length(poss_eps_sizes)
    all_ball_vars(:,:,ee) = arrayfun(@(x) x.eps_ball_var(ee),all_Cdata);
end

all_ball_alphas = bsxfun(@rdivide,all_psth_vars,all_ball_vars);
all_spline_alphas = all_psth_vars./all_spline_vars;
all_spline_alphas_noLOO = all_psth_vars./all_spline_vars_noLOO;

all_mod_psth_vars = arrayfun(@(x) mean(x.mod_psth_vars),all_Cdata);
all_mod_tot_vars = arrayfun(@(x) mean(x.mod_tot_vars),all_Cdata);
all_mod_alphas = all_mod_psth_vars./all_mod_tot_vars;

all_psth_FF = arrayfun(@(x) x.psth_FF,all_Cdata);
all_spline_FF = arrayfun(@(x) x.spline_FF,all_Cdata);

all_mod_xvLLs = arrayfun(@(x) x.bestGQM.xvLLimp,all_Cdata(:,1));
all_null_xvLLs = arrayfun(@(x) x.nullMod.xvLLimp,all_Cdata(:,1));
all_mod_xvLLimps = all_mod_xvLLs - all_null_xvLLs;

all_SU_Lratio = arrayfun(@(x) x.unit_data.SU_Lratio,all_Cdata(:,1));
all_SU_isodist = arrayfun(@(x) x.unit_data.SU_isodist,all_Cdata(:,1));
all_SU_dprime = arrayfun(@(x) x.unit_data.SU_dprime,all_Cdata(:,1));
all_SU_rate_stability = arrayfun(@(x) x.unit_data.rate_stability_cv,all_Cdata(:,1));

RF_ecc = arrayfun(@(x) x.tune_props.RF_ecc,all_Cdata(:,1));
RF_width = 2*arrayfun(@(x) x.tune_props.RF_sigma,all_Cdata(:,1));
% RF_PSF = arrayfun(@(x) x.tune_props.RF_FSF,all_Cdata(:,1));
RF_PSF = arrayfun(@(x) x.tune_props.RF_gSF,all_Cdata(:,1));
RF_PRM = arrayfun(@(x) x.tune_props.PRM,all_Cdata(:,1));

all_monkey = {all_Cdata(:,1).monkey};
all_CID = [all_Cdata(:,1).cell_ID];

pair_IDs = cat(1,all_Pdata(:,1).cell_IDs);
pair_matches = nan(size(pair_IDs));
for ii = 1:size(pair_IDs,1)
    pair_matches(ii,1) = find(all_CID == pair_IDs(ii,1));
    pair_matches(ii,2) = find(all_CID == pair_IDs(ii,2));
end


min_nTrials = 25;
min_avgRate = 5;
min_xvLL = 0;
uset = find(all_ntrials >= min_nTrials & all_avgrates >= min_avgRate & all_mod_xvLLimps > min_xvLL);
upairs = find(all(ismember(pair_IDs,all_CID(uset)),2) & pair_IDs(:,1) ~= pair_IDs(:,2));
upairs_acorr = find(all(ismember(pair_IDs,all_CID(uset)),2) & pair_IDs(:,1) == pair_IDs(:,2));

fprintf('Using %d SUs, %d pairs\n',length(uset),length(upairs));

%% compare model-predicted and direct estimates of alpha
close all

mSize = 10;
dt_ind = 1;

f1 = figure(); hold on
jbe_units = uset(strcmp(all_monkey(uset),'jbe'));
lem_units = uset(strcmp(all_monkey(uset),'lem'));
% plot(all_mod_alphas(lem_units),all_spline_alphas(lem_units),'o');
% plot(all_mod_alphas(jbe_units),all_spline_alphas(jbe_units),'ro');

% plot(all_mod_alphas(lem_units,dt_ind),all_ball_alphas(lem_units,1,dt_ind),'.','markersize',mSize);
% plot(all_mod_alphas(jbe_units,dt_ind),all_ball_alphas(jbe_units,1,dt_ind),'r.','markersize',mSize);
plot(all_mod_alphas(lem_units,dt_ind),all_spline_alphas(lem_units,dt_ind),'.','markersize',mSize);
plot(all_mod_alphas(jbe_units,dt_ind),all_spline_alphas(jbe_units,dt_ind),'r.','markersize',mSize);
line([0 1],[0 1],'color','k');
xlabel('Model-predicted alpha');
ylabel('Direct estimate alpha');
legend('LEM','JBE','Location','Southeast');

% [a,b] = corr(all_mod_alphas(uset,dt_ind),all_ball_alphas(uset,1,dt_ind),'type','pearson');
[a,b] = corr(all_mod_alphas(uset,dt_ind),all_spline_alphas(uset,dt_ind),'type','pearson');
title(sprintf('corr: %.3f\n',a));

% fig_width = 4; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Model_vs_direct_alpha.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);

%% DIRECT ESTIMATES OF ALPHA VS RF PROPERTIES
close all

mSize = 10;

f1 = figure();
subplot(2,2,1)
% plot(RF_ecc(uset),all_spline_alphas(uset),'o')
plot(RF_ecc(uset),all_ball_alphas(uset,1),'.','markersize',mSize)
[a,b] = corr(RF_ecc(uset),all_ball_alphas(uset,1),'type','spearman');
title(sprintf('ECC corr; %.3f, p %.2g\n',a,b));
xlabel('Eccentricity (deg)');
xlim([0 5]);

subplot(2,2,2)
% plot(RF_width(uset),all_spline_alphas(uset),'o')
plot(RF_width(uset),all_ball_alphas(uset,1),'.','markersize',mSize)
set(gca,'xscale','log'); xlim([0.075 1.5])
[a,b] = corr(RF_width(uset),all_ball_alphas(uset,1),'type','spearman');
title(sprintf('Width corr; %.3f, p %.2g\n',a,b));
xlabel('RF width (deg)');

subplot(2,2,3)
% plot(RF_PSF(uset),all_spline_alphas(uset),'o')
plot(RF_PSF(uset),all_ball_alphas(uset,1),'.','markersize',mSize)
[a,b] = corr(RF_PSF(uset),all_ball_alphas(uset,1),'type','spearman');
title(sprintf('SF corr; %.3f, p %.2g\n',a,b));
xlabel('Preferred SF (cyc/deg)');

subplot(2,2,4)
% plot(RF_PRM(uset),all_spline_alphas(uset),'o')
plot(RF_PRM(uset),all_ball_alphas(uset,1),'.','markersize',mSize)
[a,b] = corr(RF_PRM(uset),all_ball_alphas(uset,1),'type','spearman');
title(sprintf('PRM corr; %.3f, p %.2g\n',a,b));
xlabel('PRM');

% fig_width = 8; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Direct_alpha_vs_RF.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);

%% COMPARE DIRECT FF ESTIMATES
close all

mSize = 10;
dt_ind = 3;

f1 = figure();
plot(all_psth_FF(uset,dt_ind),all_spline_FF(uset,dt_ind),'.','markersize',mSize)
line([0 2],[0 2],'color','k');
line([0 2],[1 1],'color','k','linestyle','--');
line([1 1],[0 2],'color','k','linestyle','--');
xlabel('PSTH-based FF');
ylabel('EP-corrected FF');

% SNR = all_ball_vars(uset,1)./(all_avgrates(uset)*EP_params.base_dt);
% FF_diff = all_psth_FF(uset) - all_spline_FF(uset);
% f2 = figure();
% plot(SNR,FF_diff,'.','markersize',mSize);
% xlabel('SNR');
% ylabel('FF difference');

% fig_width = 4; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Direct_FF_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);
% 
% fig_width = 4; rel_height = 1;
% figufy(f2);
% fname = [fig_dir 'Direct_FF_SNR_compare.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);

%% direct acorr estimation
close all

% all_EP_acovs = cat(1,cell2mat(arrayfun(@(x) mean(x.EP_xcovar,1),all_Pdata(upairs_acorr),'uniformoutput',0)));
all_EP_acovs = cat(1,cell2mat(arrayfun(@(x) mean(x.EP_xcovar_LOO,1),all_Pdata(upairs_acorr),'uniformoutput',0)));
all_psth_acovs = cat(1,cell2mat(arrayfun(@(x) mean(x.pair_xcovar,1),all_Pdata(upairs_acorr),'uniformoutput',0)));
all_tot_acovs = cat(1,cell2mat(arrayfun(@(x) mean(x.tot_xcovar,1),all_Pdata(upairs_acorr),'uniformoutput',0)));
 dt = EP_params.base_dt;
all_EP_acorrs = bsxfun(@rdivide,all_EP_acovs,all_tot_acovs(:,tlags==0));
all_psth_acorrs = bsxfun(@rdivide,all_psth_acovs,all_tot_acovs(:,tlags==0));
%  all_EP_acorrs = bsxfun(@rdivide,all_EP_acovs,all_EP_acovs(:,tlags==0));
% all_psth_acorrs = bsxfun(@rdivide,all_psth_acovs,all_psth_acovs(:,tlags==0));

f1 = figure(); hold on
shadedErrorBar(tlags*dt,nanmean(all_EP_acorrs),nanstd(all_EP_acorrs)/sqrt(length(upairs_acorr)));
shadedErrorBar(tlags*dt,nanmean(all_psth_acorrs),nanstd(all_psth_acorrs)/sqrt(length(upairs_acorr)),{'color','r'});
xlim([0 0.1]);
line([0 0.1],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)');
ylabel('Autocorrelation');

% fig_width = 4; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Direct_acorr.pdf'];
% exportfig(1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);

%% xcorr analysis
close all

all_Cvars = arrayfun(@(x) mean(x.tot_var),all_Cdata);

all_tot_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.tot_xcovar,1),all_Pdata(upairs,:),'uniformoutput',0)));
% bad_pairs = find(isnan(all_tot_xcovs(:,11)));
bad_pairs = find(isnan(all_tot_xcovs(:,1)));
upairs(bad_pairs) = []; all_tot_xcovs(bad_pairs,:) = [];

all_psth_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.pair_xcovar,1),all_Pdata(upairs,:),'uniformoutput',0)));
% all_EP_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.EP_xcovar,1),all_Pdata(upairs),'uniformoutput',0)));
% all_EP_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.EP_xcovar,1),all_Pdata(upairs),'uniformoutput',0)));
all_EP_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.EP_xcovar_LOO,1),all_Pdata(upairs,:),'uniformoutput',0)));

all_norms = sqrt(all_Cvars(pair_matches(upairs,1),:).*all_Cvars(pair_matches(upairs,2),:));

all_psth_noisecorr = bsxfun(@rdivide,all_tot_xcovs - all_psth_xcovs,all_norms);
all_EP_noisecorr = bsxfun(@rdivide,all_tot_xcovs - all_EP_xcovs,all_norms);
all_psth_sigcorr = bsxfun(@rdivide,all_psth_xcovs,all_norms);
all_EP_sigcorr = bsxfun(@rdivide,all_EP_xcovs,all_norms);

xl1 = [-0.3 0.3]; 
xl2 = [-0.1 0.1];
mSize = 8;

f1 = figure();
subplot(2,2,1);hold on
plot(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11),'.','markersize',mSize)
r1 = robustfit(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11));
xx = linspace(-0.3,0.3,100);
plot(xx,r1(1) + r1(2)*xx,'k')
xlim(xl1); ylim(xl1);
line(xl1,[0 0],'color','k','linestyle','--'); line([0 0],xl1,'color','k','linestyle','--');
xlabel('Signal corr');
ylabel('Noise corr');
title('PSTH-based');

subplot(2,2,2);hold on
plot(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11),'r.','markersize',mSize)
r2 = robustfit(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11));
plot(xx,r2(1) + r2(2)*xx,'g')
xlim(xl1); ylim(xl1);
line(xl1,[0 0],'color','k','linestyle','--'); line([0 0],xl1,'color','k','linestyle','--');
xlabel('Signal corr');
ylabel('Noise corr');
title('EP-corrected');

subplot(2,2,3);hold on
plot(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11),'.','markersize',mSize)
r1 = robustfit(all_psth_sigcorr(:,11),all_psth_noisecorr(:,11));
xx = linspace(-0.2,0.2,100);
plot(xx,r1(1) + r1(2)*xx,'k')
xlim(xl2); ylim(xl2);
line(xl1,[0 0],'color','k','linestyle','--'); line([0 0],xl1,'color','k','linestyle','--');
xlabel('Signal corr');
ylabel('Noise corr');
title('PSTH-based');

subplot(2,2,4);hold on
plot(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11),'r.','markersize',mSize)
r2 = robustfit(all_EP_sigcorr(:,11),all_EP_noisecorr(:,11));
plot(xx,r2(1) + r2(2)*xx,'g')
xlim(xl2); ylim(xl2);
line(xl1,[0 0],'color','k','linestyle','--'); line([0 0],xl1,'color','k','linestyle','--');
xlabel('Signal corr');
ylabel('Noise corr');
title('EP-corrected');


dt = EP_params.base_dt;
 sig_yl = [-0.05 0.075];
 EP_base_sigcorrs = all_EP_sigcorr(:,tlags==0);
 %  negcorr_set = find(EP_base_sigcorrs <= prctile(EP_base_sigcorrs,25));
 %  poscorr_set = find(EP_base_sigcorrs >= prctile(EP_base_sigcorrs,75));
 cthresh = 0.025;
 negcorr_set = find(EP_base_sigcorrs <= -cthresh);
 poscorr_set = find(EP_base_sigcorrs >= cthresh);
 f2 = figure();
 subplot(2,2,1);
 shadedErrorBar(tlags*dt,nanmean(all_psth_sigcorr(poscorr_set,:)),nanstd(all_psth_sigcorr(poscorr_set,:))/sqrt(length(poscorr_set)),{'color','r'});
 hold on
  shadedErrorBar(tlags*dt,nanmean(all_psth_sigcorr(negcorr_set,:)),nanstd(all_psth_sigcorr(negcorr_set,:))/sqrt(length(negcorr_set)),{'color','b'});
 ylim(sig_yl);
 line([-0.1 0.1],[0 0],'color','k');
 title('PSTH signal correlation');
 xlabel('Time (s)');
 ylabel('Correlation');
 subplot(2,2,2);
 shadedErrorBar(tlags*dt,nanmean(all_EP_sigcorr(poscorr_set,:)),nanstd(all_EP_sigcorr(poscorr_set,:))/sqrt(length(poscorr_set)),{'color','r'});
 hold on
  shadedErrorBar(tlags*dt,nanmean(all_EP_sigcorr(negcorr_set,:)),nanstd(all_EP_sigcorr(negcorr_set,:))/sqrt(length(negcorr_set)),{'color','b'});
 ylim(sig_yl);
 line([-0.1 0.1],[0 0],'color','k');
 title('EP-corrected signal correlation');
  xlabel('Time (s)');
 ylabel('Correlation');
 subplot(2,2,3);
 shadedErrorBar(tlags*dt,nanmean(all_psth_noisecorr(poscorr_set,:)),nanstd(all_psth_noisecorr(poscorr_set,:))/sqrt(length(poscorr_set)),{'color','r'});
 hold on
  shadedErrorBar(tlags*dt,nanmean(all_psth_noisecorr(negcorr_set,:)),nanstd(all_psth_noisecorr(negcorr_set,:))/sqrt(length(negcorr_set)),{'color','b'});
 ylim(sig_yl);
 line([-0.1 0.1],[0 0],'color','k');
 title('PSTH noise correlation');
 xlabel('Time (s)');
 ylabel('Correlation');
 subplot(2,2,4);
 shadedErrorBar(tlags*dt,nanmean(all_EP_noisecorr(poscorr_set,:)),nanstd(all_EP_noisecorr(poscorr_set,:))/sqrt(length(poscorr_set)),{'color','r'});
 hold on
  shadedErrorBar(tlags*dt,nanmean(all_EP_noisecorr(negcorr_set,:)),nanstd(all_EP_noisecorr(negcorr_set,:))/sqrt(length(negcorr_set)),{'color','b'});
 ylim(sig_yl);
 line([-0.1 0.1],[0 0],'color','k');
 title('EP-corrected noise correlation');
 xlabel('Time (s)');
 ylabel('Correlation');
  
  %  utlags = find(abs(tlags) <= 3);
%  psth_slope = nan(length(upairs),1);
%  ep_slope = nan(length(upairs),1);
%  for ii = 1:length(upairs)
%      psth_slope(ii) = regress(all_psth_noisecorr(ii,utlags)',all_psth_sigcorr(ii,utlags)');
%      ep_slope(ii) = regress(all_EP_noisecorr(ii,utlags)',all_EP_sigcorr(ii,utlags)');
%  end
 
fig_width = 8; rel_height = 1;
figufy(f1);
fname = [fig_dir 'Xcorr_scatter.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

fig_width = 8; rel_height = 1;
figufy(f2);
fname = [fig_dir 'Xcorr_functions.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
close all

pair_id = 161;
Expt_num = all_Pdata(pair_id).Expt_num;
bar_ori = all_Pdata(pair_id).bar_ori;
dname = sprintf('~/Analysis/bruce/M%.3d/models/corrected_models_comp_ori%d',Expt_num,bar_ori)
temp = load(dname,'modFitParams');
sp_dx = temp.modFitParams.sp_dx;
mod_dt = temp.modFitParams.dt;
dt = EP_params.base_dt;

f1 = figure(); 
subplot(2,1,1);
hold on
plot(tlags*dt,mean(all_Pdata(pair_id).tot_xcovar,1),'k','linewidth',2);
plot(tlags*dt,mean(all_Pdata(pair_id).pair_xcovar,1),'r','linewidth',2);
plot(tlags*dt,mean(all_Pdata(pair_id).EP_xcovar,1),'b','linewidth',2);
xlabel('Time (s)');
ylabel('Correlation');
subplot(2,1,2);
hold on
plot(tlags*dt,mean(all_Pdata(pair_id).mod_tot_covar,1),'k','linewidth',2);
plot(tlags*dt,mean(all_Pdata(pair_id).mod_psth_covar,1),'r','linewidth',2);
xlabel('Time (s)');
ylabel('Correlation');

xr = [-0.4 0.4];

c1 = find([all_Cdata(:).cell_ID] == all_Pdata(pair_id).cell_IDs(1));
c2 = find([all_Cdata(:).cell_ID] == all_Pdata(pair_id).cell_IDs(2));
[f2,c1_dims,mod_filts1,mod_signs1] = plot_mod_filters(all_Cdata(c1).bestGQM,sp_dx,mod_dt);
[f3,c2_dims,mod_filts2,mod_signs2] = plot_mod_filters(all_Cdata(c2).bestGQM,sp_dx,mod_dt);
figure(f2);
ch = get(f2,'children');
for ii = 1:length(ch)
    if strcmp(get(ch(ii),'type'),'axes')
    xlim(ch(ii),xr);
    end
end
c1_n_wins = length(ch);

ch = get(f3,'children');
for ii = 1:length(ch)
    if strcmp(get(ch(ii),'type'),'axes')
    xlim(ch(ii),xr);
    end
end
c2_n_wins = length(ch);

eq1 = find(mod_signs1(2:end) == 1) + 1;
eq2 = find(mod_signs2(2:end) == 1) + 1;
nPix = size(mod_filts1,2); flen = size(mod_filts1,1);
xax = (1:nPix)*sp_dx; xax = xax - mean(xax);
tax = (0:(flen-1))*dt + dt/2; tax = tax*1e3;
f4 = figure();
subplot(2,1,1);
imagesc(xax,tax,sqrt(squeeze(sum(mod_filts1(:,:,eq1).^2,3))));
set(gca,'ydir','normal'); xlim(xr);
xlabel('Rel position (deg)');
ylabel('Time lag (ms)');
subplot(2,1,2);
imagesc(xax,tax,sqrt(squeeze(sum(mod_filts2(:,:,eq2).^2,3))));
set(gca,'ydir','normal'); xlim(xr);
xlabel('Rel position (deg)');
ylabel('Time lag (ms)');

fig_width = 4; rel_height = 1.6;
figufy(f1);
fname = [fig_dir sprintf('Xcorr_examp%d.pdf',pair_id)];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 3*c1_dims(2); rel_height = c1_dims(1)/c1_dims(2)*0.9;
% figufy(f2);
% fname = [fig_dir sprintf('Xcorr_examp%d_mod1.pdf',pair_id)];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f2);
% 
% fig_width = 3*c2_dims(2); rel_height = c2_dims(1)/c2_dims(1)*0.9;
% figufy(f3);
% fname = [fig_dir sprintf('Xcorr_examp%d_mod2.pdf',pair_id)];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f3);

fig_width = 4; rel_height = 2;
figufy(f4);
fname = [fig_dir sprintf('Xcorr_examp%d_efp.pdf',pair_id)];
exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% load in model-based calculations

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
% Expt_list = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};
% expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
expt_mname = repmat({'lem'},1,length(Expt_list));
expt_rnum = ones(length(Expt_list),2);

Expt_list = cat(2,Expt_list,{'G085','G086','G087','G088','G089','G091','G093','G095'});
expt_oris = cat(1,expt_oris,[0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan]);
expt_mname = cat(2,expt_mname,repmat({'jbe'},1,8));
expt_rnum = cat(1,expt_rnum,ones(8,2));

Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012','M013'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 160; 0 0; 100 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe','jbe'});
expt_rnum = cat(1,expt_rnum,[1 1; 1 1; 1 1; 1 1; 1 2; 1 2; 1 1]);

base_sname = 'model_variability_compact';
base_gname = 'grating_sim';

Mcnt = 1;
all_Mdata = [];
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

            gname = [data_dir base_gname sprintf('_ori%d',bar_ori)];
            if rec_number > 1
                gname = strcat(gname,sprintf('_r%d',rec_number));
            end
            gdat = load(gname);
            
            EP_data = EP_data(targs);
            for cc = 1:length(EP_data)
                if ~isempty(EP_data(cc).ModData.unit_data)
                    fprintf('Cell %d/%d\n',cc,length(EP_data));
                    
                    EP_data(cc).monkey = monk_name;
                    EP_data(cc).Expt_num = str2num(Expt_name(2:end));
                    EP_data(cc).bar_ori = bar_ori;
                    EP_data(cc).rec_number = rec_number;
                    EP_data(cc).cell_ID = Ccnt;
                    EP_data(cc).ov_EP_xcov = ov_EP_data.EP_xcov;
                    EP_data(cc).ov_EP_lags = ov_EP_data.EP_lags;
                    
                    EP_data(cc).grate_data = gdat.grate_Cdata(cc);
                    EP_data(cc).grate_ubins = gdat.poss_ubins;
                    
                    all_Mdata = cat(1,all_Mdata,EP_data(cc));
                    Mcnt = Mcnt + 1;                    
                end
            end
            
        end
    end
end

%
SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,all_Mdata);
Expt_numbers = [all_Mdata(:).Expt_num]';
Rec_numbers = [all_Mdata(:).rec_number]';
to_eliminate = [];
for ii = 1:length(all_Mdata)
    curset = find(SU_numbers == SU_numbers(ii) & Expt_numbers == Expt_numbers(ii) & Rec_numbers == Rec_numbers(ii));
    if length(curset) > 1
        cur_xvLLs = arrayfun(@(x) x.ModData.bestGQM.xvLLimp,all_Mdata(curset));
        avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_Mdata(curset));
        xvLL_rate = cur_xvLLs.*avg_rates;
        
        [~,best_ind] = max(cur_xvLLs);
%         [~,best_ind] = max(avg_rates);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
    end
end
to_eliminate = unique(to_eliminate);
fprintf('Eliminating %d/%d duplicate SUs\n',length(to_eliminate),length(all_Mdata));
all_Mdata(to_eliminate) = [];

%
SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,all_Mdata);
Expt_numbers = [all_Mdata(:).Expt_num]';
Rec_numbers = [all_Mdata(:).rec_number]';
to_eliminate = [];
for ii = 1:size(dup_SUs,1)
   cur_unit_1 = find(Expt_numbers == dup_SUs(ii,1) & Rec_numbers == 2 & SU_numbers == dup_SUs(ii,2));
   cur_unit_2 = find(Expt_numbers == dup_SUs(ii,1) & Rec_numbers == 1 & SU_numbers == dup_SUs(ii,3));
   curset = [cur_unit_1 cur_unit_2];
   if length(curset) == 2
        cur_xvLLs = arrayfun(@(x) x.ModData.bestGQM.xvLLimp,all_Mdata(curset));
        avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_Mdata(curset));
        xvLL_rate = cur_xvLLs.*avg_rates;
%         [~,best_ind] = max(avg_rates);
        [~,best_ind] = max(xvLL_rate);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
   end
end

double_CIDs = [all_Mdata(to_eliminate).cell_ID];
fprintf('Eliminating %d/%d duplicate SUs (multiple recs)\n',length(to_eliminate),length(all_Mdata));
elim_CIDs = [all_Mdata(to_eliminate).cell_ID];
all_Mdata(to_eliminate) = [];

%% select cells for analysis
all_avgrates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_Mdata);
all_monkey = {all_Mdata(:).monkey};
all_CID = [all_Mdata(:).cell_ID];

RF_ecc = arrayfun(@(x) x.ModData.tune_props.RF_ecc,all_Mdata);
RF_width = 2*arrayfun(@(x) x.ModData.tune_props.RF_sigma,all_Mdata);
% RF_PSF = arrayfun(@(x) x.ModData.tune_props.RF_FSF,all_Mdata);
RF_PSF = arrayfun(@(x) x.ModData.tune_props.RF_gSF,all_Mdata);
RF_PRM = arrayfun(@(x) x.ModData.tune_props.PRM,all_Mdata);

actual_EP_SDs = arrayfun(@(x) x.poss_SDs(end),all_Mdata);

mod_alpha_funs = cat(1,all_Mdata.alpha_funs);

min_avgRate = 5;
% uset = find(all_avgrates >= min_avgRate & RF_ecc > 1);
MD_uset = find(all_avgrates >= min_avgRate);

%% plot model-predicted alphas vs RF properties
close all
poss_SDs = all_Mdata(1).poss_SDs;

SD_ind = 3; %use this value for EP SD
if SD_ind == length(poss_SDs)
    fprintf('Evaluated at native EP SD\n');
else
    fprintf('Evaluated at EP SD of %.2f\n',poss_SDs(SD_ind));
end

mSize = 10;

f1 = figure();
subplot(2,2,1)
plot(RF_ecc(MD_uset),mod_alpha_funs(MD_uset,SD_ind),'.','markersize',mSize)
xlim([0 5]);
[a,b] = corr(RF_ecc(MD_uset),mod_alpha_funs(MD_uset,SD_ind),'type','spearman');
title(sprintf('corr; %.3f, p %.2g\n',a,b));
xlabel('Eccentricity (deg)');

subplot(2,2,2)
plot(RF_width(MD_uset),mod_alpha_funs(MD_uset,SD_ind),'.','markersize',mSize)
set(gca,'xscale','log'); xlim([0.07 1.5])
[a,b] = corr(RF_width(MD_uset),mod_alpha_funs(MD_uset,SD_ind),'type','spearman');
title(sprintf('corr; %.3f, p %.2g\n',a,b));
xlabel('RF width (deg)');

subplot(2,2,3)
plot(RF_PSF(MD_uset),mod_alpha_funs(MD_uset,SD_ind),'.','markersize',mSize)
[a,b] = corr(RF_PSF(MD_uset),mod_alpha_funs(MD_uset,SD_ind),'type','spearman');
title(sprintf('corr; %.3f, p %.2g\n',a,b));
xlabel('Preferred SF (cyc/deg)');

subplot(2,2,4)
plot(RF_PRM(MD_uset),mod_alpha_funs(MD_uset,SD_ind),'.','markersize',mSize)
[a,b] = corr(RF_PRM(MD_uset),mod_alpha_funs(MD_uset,SD_ind),'type','spearman');
title(sprintf('corr; %.3f, p %.2g\n',a,b));
xlabel('PRM');


f2 = figure();
plot(RF_ecc(MD_uset),RF_width(MD_uset),'.','markersize',mSize)
set(gca,'yscale','log'); ylim([0.07 1.5])
xlim([0 5]);
[a,b] = corr(RF_ecc(MD_uset),RF_width(MD_uset),'type','spearman');
title(sprintf('corr; %.3f, p %.2g\n',a,b));
xlabel('Eccentricity (deg)');
ylabel('RF width (deg)');

fig_width = 8; rel_height = 1;
figufy(f1);
fname = [fig_dir 'Model_alpha_vs_RF.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

fig_width = 4; rel_height = 1;
figufy(f2);
fname = [fig_dir 'RF_width_vs_ecc.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% Average eye position acorr function
close all

ov_EP_xcov = cat(2,all_Mdata(MD_uset).ov_EP_xcov);
EP_lags = all_Mdata(1).ov_EP_lags;

f1 = figure;
shadedErrorBar(EP_lags*.01,nanmean(ov_EP_xcov,2),nanstd(ov_EP_xcov,[],2));
xlabel('Time lag (s)');
ylabel('Correlation');
xlim([0 0.5]);

fig_width = 4; rel_height = 1;
figufy(f1);
fname = [fig_dir 'EP_acorr_fun.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% plot Poisson Fano-factors for different time binnings
poss_ubins = all_Mdata(1).poss_ubins;
poss_SDs = all_Mdata(1).poss_SDs(1:end-1);

ep_FF_funs = cat(3,all_Mdata.simrpt_FF);
ep_PSTH_funs = cat(3,all_Mdata.bin_PSTH_vars);
ep_tvar_funs = cat(3,all_Mdata.bin_tot_vars);
% ep_tavg_funs = cat(3,all_Mdata.bin_tot_avgs);
ep_tavg_funs = bsxfun(@times,all_avgrates',poss_ubins');
ep_tavg_funs = reshape(ep_tavg_funs,length(poss_ubins),1,[]);
ep_tavg_funs = repmat(ep_tavg_funs,[1 length(poss_SDs)+1 1]);

ep_alpha_funs = ep_PSTH_funs./ep_tvar_funs;
ep_SNR_funs = ep_tvar_funs./ep_tavg_funs;

rate_covs = arrayfun(@(x) x.base_vars(end),all_Mdata);
rate_FF = rate_covs./all_avgrates;

f1 = figure();
hold on
% plot(poss_ubins*.01,squeeze(ep_FF_funs(:,end,MD_uset)),'k','linewidth',0.5);
errorbar(poss_ubins*0.01,squeeze(nanmean(ep_FF_funs(:,end,MD_uset),3)),squeeze(nanstd(ep_FF_funs(:,end,MD_uset),[],3))/sqrt(length(MD_uset)),'color','r','linewidth',3)
xlim([0 1]);

SNR = squeeze(ep_SNR_funs(1,end,MD_uset));
hSNR_set = MD_uset(SNR > prctile(SNR,75));
lSNR_set = MD_uset(SNR < prctile(SNR,25));
errorbar(poss_ubins*0.01,squeeze(nanmean(ep_FF_funs(:,end,hSNR_set),3)),squeeze(nanstd(ep_FF_funs(:,end,hSNR_set),[],3))/sqrt(length(hSNR_set)),'color','b','linewidth',2)
errorbar(poss_ubins*0.01,squeeze(nanmean(ep_FF_funs(:,end,lSNR_set),3)),squeeze(nanstd(ep_FF_funs(:,end,lSNR_set),[],3))/sqrt(length(lSNR_set)),'color','k','linewidth',2)

xlabel('Bin width (s)');
ylabel('Fano-factor');

fig_width = 4; rel_height = 1;
figufy(f1);
fname = [fig_dir 'Model_FF_binning.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% dependence of model-predicted alpha on EP SD
close all

poss_SDs = all_Mdata(1).poss_SDs(1:end-1);
uSD_ind = 2;
mid_values = mod_alpha_funs(MD_uset,uSD_ind);
strong_set = MD_uset(mid_values < prctile(mid_values,25));
weak_set = MD_uset(mid_values > prctile(mid_values,75));

f1 = figure(); 
hold on
errorbar(poss_SDs,nanmean(mod_alpha_funs(MD_uset,1:end-1)),nanstd(mod_alpha_funs(MD_uset,1:end-1))/sqrt(length(MD_uset)));
errorbar(poss_SDs,nanmean(mod_alpha_funs(strong_set,1:end-1)),nanstd(mod_alpha_funs(strong_set,1:end-1))/sqrt(length(strong_set)),'color','r');
errorbar(poss_SDs,nanmean(mod_alpha_funs(weak_set,1:end-1)),nanstd(mod_alpha_funs(weak_set,1:end-1))/sqrt(length(weak_set)),'color','k');
xlim([0 0.2])
xlabel('Eye position SD (deg)');
ylabel('Alpha');

f3 = figure(); 
hold on
plot(poss_SDs,mod_alpha_funs(MD_uset,1:end-1),'r','linewidth',0.5);
errorbar(poss_SDs,nanmean(mod_alpha_funs(MD_uset,1:end-1)),nanstd(mod_alpha_funs(MD_uset,1:end-1)),'k','linewidth',3);
xlim([0 0.2])
xlabel('Eye position SD (deg)');
ylabel('Alpha');
ylim([0 1]);

f2 = figure(); 
subplot(2,1,1);
hist(actual_EP_SDs(MD_uset),20);
xlim([0 0.2]);
xlabel('Eye position SD (deg)');
ylabel('Number of units');

% avg_sp_fit = spline(poss_SDs,nanmean(mod_alpha_funs(MD_uset,1:end-1)),actual_EP_SDs(MD_uset));
% strong_sp_fit = spline(poss_SDs,nanmean(mod_alpha_funs(strong_set,1:end-1)),actual_EP_SDs(MD_uset));
% weak_sp_fit = spline(poss_SDs,nanmean(mod_alpha_funs(weak_set,1:end-1)),actual_EP_SDs(MD_uset));
% cur_bin_edges = linspace(0,1,50);
% subplot(2,1,2); hold on
% stairs(cur_bin_edges,histc(avg_sp_fit,cur_bin_edges)/length(MD_uset),'b');
% stairs(cur_bin_edges,histc(strong_sp_fit,cur_bin_edges)/length(MD_uset),'r');
% stairs(cur_bin_edges,histc(weak_sp_fit,cur_bin_edges)/length(MD_uset),'k');
% xlim([0 1]);
% xlabel('Alpha');
% ylabel('Relative frequency');

alpha_CVs = nan(length(MD_uset),1);
for  ii = 1:length(MD_uset)
    sp_vals = spline(poss_SDs,mod_alpha_funs(MD_uset(ii),1:end-1),actual_EP_SDs(MD_uset));
    alpha_CVs(ii) = nanstd(sp_vals)/nanmean(sp_vals);
end
subplot(2,1,2);
plot(mod_alpha_funs(MD_uset,poss_SDs==0.1),alpha_CVs,'.','markersize',10);
xlim([0 1]);
xlabel('Alpha');
ylabel('Alpha CV');


% fig_width = 4; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Alpha_vs_EPSD.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);
% 
fig_width = 4; rel_height = 2;
figufy(f2);
fname = [fig_dir 'Alpha_vs_EPSD_dists.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);


fig_width = 4; rel_height = 1;
figufy(f3);
fname = [fig_dir 'Alpha_vs_EPSD_shade.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%% Grating simulation fano factors
close all

poss_grate_sf = [1 2 4];
poss_grate_tf = [2 4 8];

grate_ubins = all_Mdata(1).grate_ubins;
all_grate_alphas = nan(length(MD_uset),length(grate_ubins),length(poss_grate_sf),length(poss_grate_tf));
all_grate_FFs = nan(length(MD_uset),length(grate_ubins),length(poss_grate_sf),length(poss_grate_tf));
all_grate_PD = nan(length(MD_uset),length(poss_grate_sf),length(poss_grate_tf));
for sf_ind = 1:length(poss_grate_sf)
    for tf_ind = 1:length(poss_grate_tf)
        for ii = 1:length(MD_uset)
            d1_var = all_Mdata(MD_uset(ii)).grate_data.tot_vars_NS(sf_ind,tf_ind,1);
            d2_var = all_Mdata(MD_uset(ii)).grate_data.tot_vars_NS_R(sf_ind,tf_ind,1);
            if d1_var > d2_var
                PSTH_vars = all_Mdata(MD_uset(ii)).grate_data.PSTH_vars_NS(sf_ind,tf_ind,:);
                tot_vars = all_Mdata(MD_uset(ii)).grate_data.tot_vars_NS(sf_ind,tf_ind,:);
                all_grate_FFs(ii,:,sf_ind,tf_ind) = all_Mdata(MD_uset(ii)).grate_data.FF_ests_NS(sf_ind,tf_ind,:);
                all_grate_PD(ii,sf_ind,tf_ind) = 1;
            else
                PSTH_vars = all_Mdata(MD_uset(ii)).grate_data.PSTH_vars_NS_R(sf_ind,tf_ind,:);
                tot_vars = all_Mdata(MD_uset(ii)).grate_data.tot_vars_NS_R(sf_ind,tf_ind,:);
                all_grate_FFs(ii,:,sf_ind,tf_ind) = all_Mdata(MD_uset(ii)).grate_data.FF_ests_NS_R(sf_ind,tf_ind,:);
                all_grate_PD(ii,sf_ind,tf_ind) = 2;
            end
        end
    end
end

lwidths = [1 2 4];
cmap = [1 0 0; 0 0 1; 0 0 0];
f1 = figure(); hold on
for sf = 1:length(poss_grate_sf)
    for tf = 1:length(poss_grate_tf)
        plot(grate_ubins*.01,squeeze(nanmean(all_grate_FFs(:,:,sf,tf))),'linewidth',lwidths(sf),'color',cmap(tf,:));
    end
end
xlabel('Time binning (s)');
ylabel('Fano factor');

fig_width = 4; rel_height = 1;
figufy(f1);
fname = [fig_dir 'Grating_FF.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%% 
SU_numbers = arrayfun(@(x) x.unit_data.SU_number,all_Cdata(uset));
Expt_numbers = [all_Cdata(uset).Expt_num]';
Rec_numbers = [all_Cdata(uset).rec_number]';
Bar_oris = [all_Cdata(uset).bar_ori]';

MD_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,all_Mdata(MD_uset));
MD_Expt_numbers = [all_Mdata(MD_uset).Expt_num]';
MD_Rec_numbers = [all_Mdata(MD_uset).rec_number]';
MD_Bar_oris = [all_Mdata(MD_uset).bar_ori]';

MD_match = nan(length(MD_uset),1);
for ii = 1:length(MD_uset)
    cur_match = find(SU_numbers == MD_SU_numbers(ii) & Expt_numbers == MD_Expt_numbers(ii) ...
        & Rec_numbers == MD_Rec_numbers(ii) & Bar_oris == MD_Bar_oris(ii));
    if ~isempty(cur_match)
       MD_match(ii) = cur_match; 
    end
end

matched_MUnits = find(~isnan(MD_match));



