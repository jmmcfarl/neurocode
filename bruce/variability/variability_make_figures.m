clear all
close all
rmpath('~/James_scripts/bruce/bruce_code/');

Expt_list = {};
expt_oris = [];
expt_mname = {};
expt_rnum = [];

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
expt_mname = repmat({'lem'},1,length(Expt_list));
expt_rnum = ones(length(Expt_list),2);

Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012','M013','M014','M320'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 160; 0 0; 100 nan; 40 nan; 100 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe','jbe','jbe','lem'});
expt_rnum = cat(1,expt_rnum,[1 1; 1 1; 1 1; 1 1; 1 2; 1 2; 1 1; 1 1; 1 1]);

fig_dir = '/home/james/Analysis/bruce/variability/figures/';

%% load repeat trial data
% base_rname = 'rpt_variability_compact_FIN'; %rpt-trial data base
% base_rname = 'rpt_variability_compact_FIN4_noxc'; %this has more time bins but no xxc
% base_rname = 'rpt_variability_compact_FIN5'; % has some time bins with xc
% base_rname = 'rpt_variability_compact_FIN_noextras'; % has some time bins with xc
% base_rname = 'rpt_variability_compact_FIN5'; % has some time bins with xc
base_rname = 'rpt_variability_compact_nFIN_noxc'; % has some time bins with xc
base_sname = 'sim_variability_compact_nFIN'; %sim-calc data base
% base_sname = 'sim_variability_compact_FIN2_noxc'; %newer sim-calc that has integral-based and no xc sim

cell_cnt = 1;
pair_cnt = 1;
all_cell_data = [];
all_pair_data = [];
all_expt_data = [];
% for Elist_cnt = 1:18 %loop over all experiments in the list
for Elist_cnt = 1:length(Expt_list) %loop over all experiments in the list
    Expt_name = Expt_list{Elist_cnt};
    monk_name = expt_mname{Elist_cnt};
    for bori_cnt = 1:2 %loop over both oris/rec nums for this expt
        bar_ori = expt_oris(Elist_cnt,bori_cnt);
        rec_number = expt_rnum(Elist_cnt,bori_cnt);
        if ~isnan(bar_ori)
            fprintf('Loading %s on Expt %s ori %d\n',base_rname,Expt_name,bar_ori);
            data_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
            
            rname = [data_dir base_rname sprintf('_ori%d',bar_ori)];
            sname = [data_dir base_sname sprintf('_ori%d',bar_ori)];
            if rec_number > 1
                rname = strcat(rname,sprintf('_r%d',rec_number));
                sname = strcat(sname,sprintf('_r%d',rec_number));
            end
            load(rname);
            load(sname);
            
            all_expt_data = cat(1,all_expt_data,rpt_data);
            
            %loop over cells for this rec
            reverseStr = '';
            for cc = 1:size(EP_data,1)
                if ~isnan(EP_data(cc,1).ov_avg_BS) %if the cell had any usable rpt data, include it
                    msg = sprintf('Cell %d/%d',cc,size(EP_data,1));
                    fprintf([reverseStr msg]);
                    reverseStr = repmat(sprintf('\b'),1,length(msg));                   
                    
                    %store meta data
                    EP_data(cc,1).monkey = monk_name;
                    EP_data(cc,1).Expt_num = str2num(Expt_name(2:end));
                    EP_data(cc,1).bar_ori = bar_ori;
                    EP_data(cc,1).rec_number = rec_number;
                    EP_data(cc,1).cell_ID = cell_cnt; %give each cell a unique integer ID
                                     
                    %add in simulated model-calcs
                    EP_data(cc,1).sim_data = sim_data(cc);
                    EP_data(cc,1).sim_params = sim_params;
                    
                    %add to cell list
                    all_cell_data = cat(1,all_cell_data,EP_data(cc,:));
                    cell_cnt = cell_cnt + 1; %oincrement cell cnter
                end
            end
            fprintf('\n');
            if EP_params.do_xcorrs %if we have cell-pair data
                cur_ucells = find(~isnan([EP_data(:,1).ov_avg_BS])); %find all usable cells in this expt
                for cc = 1:size(EP_pairs,1) %loop over unique pairs
                    if all(ismember(EP_pairs(cc,1).ids,cur_ucells)) %if both members were usable
                        
                        %store metadata
                        EP_pairs(cc,1).monkey = monk_name;
                        EP_pairs(cc,1).Expt_num = str2num(Expt_name(2:end));
                        EP_pairs(cc,1).bar_ori = bar_ori;
                        EP_pairs(cc,1).rec_number = rec_number;
                        EP_pairs(cc,1).cell_IDs = [EP_data([EP_pairs(cc,1).ids]).cell_ID]; %cell IDs for each cell in pair
                        EP_pairs(cc,1).pair_ID = pair_cnt;
                                                
                        %add simulated data
                        EP_pairs(cc,1).sim_data = sim_pairs(cc);
                        
                        %add pair to list and increment cnter
                        all_pair_data = cat(1,all_pair_data,EP_pairs(cc,:));
                        pair_cnt = pair_cnt + 1;
                    end
                end
            end
        end
    end
end

%% FOR CELLS RECORDED MULTIPLE TIMES, PICK BEST INSTANCE
SU_numbers = arrayfun(@(x) x.unit_data.SU_number,all_cell_data(:,1)); 
Expt_numbers = [all_cell_data(:,1).Expt_num]';
Rec_numbers = [all_cell_data(:,1).rec_number]';
bar_oris = [all_cell_data(:,1).bar_ori]';

to_eliminate = [];
for ii = 1:size(all_cell_data,1) %loop over SUs
    %find all times this cell was recorded with the same Rec_number
    %(different bar oris)
    curset = find(SU_numbers == SU_numbers(ii) & Expt_numbers == Expt_numbers(ii) & Rec_numbers == Rec_numbers(ii));
    if length(curset) > 1 %if there was a repeat rec
%         cur_xvLLs = arrayfun(@(x) x.bestGQM.xvLLimp,all_cell_data(curset,1)); %xval LL imp over null model
%         avg_rates = arrayfun(@(x) x.unit_data.avg_rate,all_cell_data(curset,1)); %avg rates of the cells
        cur_xvLLs = arrayfun(@(x) x.bestGQM.rptLLimp,all_cell_data(curset,1)); %xval LL imp over null model
        avg_rates = arrayfun(@(x) x.ov_avg_BS,all_cell_data(curset,1)); %avg rates of the cells
        xvLL_rate = cur_xvLLs.*avg_rates; %LL per unit time
%         [~,best_ind] = max(xvLL_rate); %find the best rec for this cell
        [~,best_ind] = max(cur_xvLLs); %find the best rec for this cell
        worst_ind = setdiff(1:length(curset),best_ind); %eliminate all others
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
    end
end
to_eliminate = unique(to_eliminate);

% FOR SAME SUS RECORDED ON MULTIPLE SESSIONS WITH DIFFERENT electrode
% depths (rec numbers). Just hard code the instances of this and manually
% pull them out
duplicate_SUs = [12 1 5; 12 3 8]; %[Expt_num r2_SU_Number r1_SU_number]
for ii = 1:size(duplicate_SUs,1)
   cur_unit_1 = find(Expt_numbers == duplicate_SUs(ii,1) & Rec_numbers == 2 & SU_numbers == duplicate_SUs(ii,2));
   cur_unit_2 = find(Expt_numbers == duplicate_SUs(ii,1) & Rec_numbers == 1 & SU_numbers == duplicate_SUs(ii,3));
   curset = [cur_unit_1 cur_unit_2];
   if length(curset) == 2
%         cur_xvLLs = arrayfun(@(x) x.bestGQM.xvLLimp,all_cell_data(curset,1));
%         avg_rates = arrayfun(@(x) x.unit_data.avg_rate,all_cell_data(curset,1));
        cur_xvLLs = arrayfun(@(x) x.bestGQM.rptLLimp,all_cell_data(curset,1));
        avg_rates = arrayfun(@(x) x.ov_avg_BS,all_cell_data(curset,1));
        xvLL_rate = cur_xvLLs.*avg_rates;
%         [~,best_ind] = max(xvLL_rate);
        [~,best_ind] = max(cur_xvLLs);
        worst_ind = setdiff(1:length(curset),best_ind);
        to_eliminate = cat(1,to_eliminate,curset(worst_ind));
   end
end

%%
fprintf('Eliminating %d/%d duplicate SUs (multiple recs)\n',length(to_eliminate),size(all_cell_data,1));
elim_CIDs = [all_cell_data(to_eliminate,1).cell_ID]; %IDs of cells being removed
all_cell_data(to_eliminate,:) = [];
if EP_params.do_xcorrs
    %remove pairs where either unit is one being eliminated
    pair_IDs = cat(1,all_pair_data(:,1).cell_IDs);
    all_pair_data(any(ismember(pair_IDs,elim_CIDs),2),:) = [];
end

%% extract SU properties and select units for analysis
n_SUs = size(all_cell_data,1); 
SU_monkey = {all_cell_data(:,1).monkey};
SU_CID = [all_cell_data(:,1).cell_ID];

direct_bin_dts = EP_params.direct_bin_dts; %range of time bins used for direct analysis
mod_bin_dts = EP_params.mod_bin_dts; %range of time bins used for model-based analysis
poss_bin_dts = EP_params.poss_bin_dts; %set of possible time bins
direct_used_dts = find(ismember(poss_bin_dts,direct_bin_dts));

%selection criteria
min_nTrials = 25; %minimum number of repeat trials
min_avgRate = 5; %minimum avg rate (Hz)
min_xvLL = -Inf; %minimum model xval LL improvement over null

SU_nTrials = arrayfun(@(x) sum(x.n_utrials),all_cell_data(:,1));
SU_avgRates = [all_cell_data(:,1).ov_avg_BS]'/poss_bin_dts(1); %compute avg rate using first time bin res
SU_mod_xvLLs = arrayfun(@(x) x.bestGQM.xvLLimp,all_cell_data(:,1));
SU_rpt_xvLL = arrayfun(@(x) x.rpt_LL - x.rpt_nullLL,all_cell_data(:,1));

% SU_uset = find(SU_nTrials >= min_nTrials & SU_avgRates >= min_avgRate & SU_mod_xvLLs > min_xvLL); %used SUs
base_uset = find(SU_nTrials >= min_nTrials & SU_avgRates >= min_avgRate); %used SUs
SU_uset = find(SU_nTrials >= min_nTrials & SU_avgRates >= min_avgRate & SU_rpt_xvLL > min_xvLL); %used SUs
if EP_params.do_xcorrs
    pair_IDs = cat(1,all_pair_data(:,1).cell_IDs);
    upairs = find(all(ismember(pair_IDs,SU_CID(SU_uset)),2) & pair_IDs(:,1) ~= pair_IDs(:,2)); %used pairs are different from each other, and both members of used SU set
    upairs_acorr = find(all(ismember(pair_IDs,SU_CID(SU_uset)),2) & pair_IDs(:,1) == pair_IDs(:,2)); %for autocorrelation analysis, these are used pairs that are self-same
    fprintf('Using %d SUs, %d pairs\n',length(SU_uset),length(upairs));
else
    fprintf('Using %d SUs\n',length(SU_uset));
end

%% extract useful SU properties
SU_monkey = {all_cell_data(SU_uset,1).monkey};
SU_exptNumber = [all_cell_data(SU_uset,1).Expt_num]';
SU_barOri = [all_cell_data(SU_uset,1).bar_ori]';
SU_CID = [all_cell_data(SU_uset,1).cell_ID];
SU_nTrials = arrayfun(@(x) sum(x.n_utrials),all_cell_data(SU_uset,1));
SU_avgRates = [all_cell_data(SU_uset,1).ov_avg_BS]'/poss_bin_dts(1); %compute avg rate using first time bin res
SU_mod_xvLLs = arrayfun(@(x) x.bestGQM.xvLLimp,all_cell_data(SU_uset,1));
SU_numbers = arrayfun(@(x) x.unit_data.SU_number,all_cell_data(SU_uset,1)); 

%RF properties
RF_ecc = arrayfun(@(x) x.tune_props.RF_ecc,all_cell_data(SU_uset,1)); %based on gabor-fit RF mean s
RF_ecc_avg = arrayfun(@(x) x.tune_props.RF_ecc_avg,all_cell_data(SU_uset,1)); %based on gaussian fit to avg filter
RF_width = 2*arrayfun(@(x) x.tune_props.RF_sigma,all_cell_data(SU_uset,1)); %based on gabor-fit RF mean s
RF_avg_width = 2*arrayfun(@(x) x.tune_props.avgRF_sigma,all_cell_data(SU_uset,1)); %based on gaussian fit to avg filter
RF_FSF = arrayfun(@(x) x.tune_props.RF_FSF,all_cell_data(SU_uset,1));
RF_gSF = arrayfun(@(x) x.tune_props.RF_gSF,all_cell_data(SU_uset,1));
RF_PRM = arrayfun(@(x) x.tune_props.PRM,all_cell_data(SU_uset,1));
RF_PRI = arrayfun(@(x) x.tune_props.PRI,all_cell_data(SU_uset,1));

EP_SDs = arrayfun(@(x) x.EP_SD,all_cell_data(SU_uset,1)); %robust SD of EP during repeats

%cluster quality measures
SU_Lratio = arrayfun(@(x) x.unit_data.SU_Lratio,all_cell_data(SU_uset,1));
SU_isodist = arrayfun(@(x) x.unit_data.SU_isodist,all_cell_data(SU_uset,1));
SU_dprime = arrayfun(@(x) x.unit_data.SU_dprime,all_cell_data(SU_uset,1));
SU_rate_stability = arrayfun(@(x) x.unit_data.rate_stability_cv,all_cell_data(SU_uset,1));

n_SUs_used = length(SU_uset);

%% get core rate variance estimates

Mod_psth_vars = arrayfun(@(x) mean(x.mod_psth_vars_noEM),all_cell_data(SU_uset,:)); %PSTH variance of model-rates (avg over recs)
Mod_tot_vars = arrayfun(@(x) mean(x.mod_tot_vars_noEM),all_cell_data(SU_uset,:)); %total rate variance of models (avg over recs)
Mod_alphas = arrayfun(@(x) mean(x.mod_alphas_noEM),all_cell_data(SU_uset,:)); %PSTH variance of model-rates (avg over recs)

SU_psth_vars = arrayfun(@(x) mean(x.pair_psth_var),all_cell_data(SU_uset,direct_used_dts)); %unbiased PSTH variance est
SU_tot_vars = arrayfun(@(x) mean(x.tot_var),all_cell_data(SU_uset,direct_used_dts)); %unbiased PSTH variance est
SU_mean_cnts = arrayfun(@(x) x.ov_avg_BS,all_cell_data(SU_uset,direct_used_dts));

%total rate variances using epsilon balls
poss_eps_sizes = EP_params.poss_eps_sizes; %possible epsilon balls
SU_ball_vars = nan(n_SUs_used,length(direct_used_dts),length(poss_eps_sizes));
SU_ball_vars_noLOO = nan(n_SUs_used,length(direct_used_dts),length(poss_eps_sizes));
for ee = 1:length(poss_eps_sizes)
    SU_ball_vars(:,:,ee) = arrayfun(@(x) x.eps_ball_var(ee),all_cell_data(SU_uset,direct_used_dts));
    SU_ball_vars_noLOO(:,:,ee) = arrayfun(@(x) x.eps_ball_var_noLOO(ee),all_cell_data(SU_uset,direct_used_dts));
end

%alpha estimates, with and without LOO
SU_ball_alphas = 1 - bsxfun(@rdivide,SU_psth_vars,SU_ball_vars);
SU_ball_alphas_noLOO = 1 - bsxfun(@rdivide,SU_psth_vars,SU_ball_vars_noLOO);

%psth and EP-based estimates of avg. noise variance
SU_psth_noisevars = SU_tot_vars - SU_psth_vars;
SU_ball_noisevars = bsxfun(@plus,-SU_ball_vars,SU_tot_vars);

%get simulation-based alphas, evaluated at rpt eye pos SD 
sim_alphas = nan(length(SU_uset),length(sim_params.poss_ubins));
sim_alphas_const = nan(length(SU_uset),length(sim_params.poss_ubins));
sim_int_alphas = nan(length(SU_uset),1);
for ii = 1:length(SU_uset)
    cur_alpha = all_cell_data(SU_uset(ii),1).sim_data.alphas;
    sim_alphas(ii,:) = interp1(sim_params.poss_SDs(1:end-1),cur_alpha(1:end-1,:),EP_SDs(ii)); %interpolate alpha onto measured EP SD during repeat trials
    if isfield(all_cell_data(SU_uset(ii),1).sim_data,'sim_int_alphas')
    cur_alpha = all_cell_data(SU_uset(ii),1).sim_data.sim_int_alphas; %this is the calculation based on the integral of the FT
    sim_int_alphas(ii) = interp1(sim_params.poss_SDs(1:end-1),cur_alpha(1:end-1),EP_SDs(ii)); %interpolate onto measured EP SD during repeat trials
    end
    if isfield(sim_params,'calc_Tconst') && sim_params.calc_Tconst
       cur_alpha = all_cell_data(SU_uset(ii),1).sim_data.alphas_const;
        sim_alphas_const(ii,:) = interp1(sim_params.poss_SDs(1:end-1),cur_alpha(1:end-1,:),EP_SDs(ii)); %interpolate alpha onto measured EP SD during repeat trials
    end
end

%% GENERAL parameter values to use for plots (dt, epsilon ball)
mSize = 10; %markersize
dt_ind = 0.01; %which dt value to use
ball_eps = 0.01; %which epsilon ball radius to use

%find indices where the desired dt resolution was computed
mod_dt_ind = find(mod_bin_dts == dt_ind);
direct_dt_ind = find(direct_bin_dts == dt_ind);
ball_ind = find(poss_eps_sizes == ball_eps);
sim_dt_ind = find(sim_params.poss_ubins*.01 == dt_ind);

%% some quantitative comparisons on different estimates of alpha
%quantify overfitting by comparing rate variance of eps-ball estimators with and without using LOO.
overfit_measure = 100*bsxfun(@rdivide,SU_ball_vars_noLOO(:,direct_dt_ind,ball_ind) - SU_ball_vars(:,direct_dt_ind,ball_ind),SU_ball_vars(:,direct_dt_ind,ball_ind));
abs_overfit_measure = 100*bsxfun(@rdivide,abs(SU_ball_vars_noLOO(:,direct_dt_ind,ball_ind) - SU_ball_vars(:,direct_dt_ind,ball_ind)),SU_ball_vars(:,direct_dt_ind,ball_ind));
% overfit_measure = SU_ball_alphas_noLOO(:,direct_dt_ind,ball_ind) - SU_ball_alphas(:,direct_dt_ind,ball_ind);
% overfit_measure = 100*bsxfun(@rdivide,SU_ball_alphas_noLOO(:,direct_dt_ind,ball_ind) - SU_ball_alphas(:,direct_dt_ind,ball_ind),SU_ball_alphas(:,direct_dt_ind,ball_ind));
fprintf('overfitting quantiles [25 50 75]: %.4f %.4f %.4f\n',prctile(overfit_measure,[25 50 75]));
fprintf('abs overfitting quantiles [25 50 75]: %.4f %.4f %.4f\n',prctile(abs_overfit_measure,[25 50 75]));

%look at simulated alpha calcs based on integrating the FT, compared to the
%normal way
sim_int_diff = sim_int_alphas - sim_alphas(:,1);
fprintf('integral mad: %.4f SD:%.4f\n',median(abs(sim_int_diff)),std(sim_int_diff));

%compare alpha simulations based on measured EPs vs assuming constant
%within-trial EP
sim_const_diff = sim_alphas_const(:,1) - sim_alphas(:,1);
fprintf('const EP mad: %.4f SD:%.4f\n',median(abs(sim_const_diff)),std(sim_const_diff));

%% compare model-predicted and direct estimates of alpha
% close all

f1 = figure(); hold on
plot(Mod_alphas(:,mod_dt_ind),SU_ball_alphas(:,direct_dt_ind,ball_ind),'.','markersize',mSize);
% plot(sim_alphas(:,sim_dt_ind),SU_ball_alphas(:,direct_dt_ind,ball_ind),'.','markersize',mSize);
line([0 1],[0 1]);
xlabel('Model alpha');
ylabel('Direct alpha');

[a,b] = corr(Mod_alphas(:,mod_dt_ind),SU_ball_alphas(:,direct_dt_ind,ball_ind),'type','pearson');
% [a,b] = corr(sim_alphas(:,sim_dt_ind),SU_ball_alphas(:,direct_dt_ind,ball_ind),'type','pearson');
title(sprintf('corr: %.3f',a));


% fig_width = 4; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Mod_vs_ball_alpha.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% plot alpha distributions
fov_cells = find(RF_ecc <= 2);
parafov_cells = find(RF_ecc > 2);

n_bins = 20;
bin_edges = linspace(0,1,n_bins + 1);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
n_fovea = histc(SU_ball_alphas(fov_cells,direct_dt_ind,ball_ind),bin_edges);
% n_fovea = n_fovea/length(fov_cells);
n_para = histc(SU_ball_alphas(parafov_cells,direct_dt_ind,ball_ind),bin_edges);
% n_para = n_para/length(parafov_cells);

% f1 = figure(); hold on
% stairs(bin_edges,n_fovea,'r');
% stairs(bin_edges,n_para,'k');

n_tot = histc(SU_ball_alphas(:,direct_dt_ind,ball_ind),bin_edges);
% n_tot = n_tot/length(SU_uset);
f2 = figure();
h = bar(gca,bin_cents,n_tot(1:end-1));
set(h,'barwidth',0.8,'faceColor','k');
xlabel('Alpha'); ylabel('Number of cells');
xlim(minmax(bin_edges));
yl = ylim();
line(median(SU_ball_alphas(:,direct_dt_ind,ball_ind))+[0 0],yl,'color','b');

% n_simple = histc(SU_ball_alphas(simple,direct_dt_ind,ball_ind),bin_edges);
% n_simple = n_simple/length(simple);
% n_complex = histc(SU_ball_alphas(complex,direct_dt_ind,ball_ind),bin_edges);
% n_complex = n_complex/length(complex);
% 
% f3 = figure(); hold on
% stairs(bin_edges,n_simple,'r');
% stairs(bin_edges,n_complex,'k');

% fig_width = 4; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Overall_alpha_dist.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% compare alphas as a fnx of time window
use_np_bins = find(poss_bin_dts <= 1);

%normalize different alpha estimates by thier values at base dt (usually
%10ms)
Mod_rel_alphas = bsxfun(@rdivide,Mod_alphas,Mod_alphas(:,mod_dt_ind));
ball_rel_alphas = bsxfun(@rdivide,SU_ball_alphas(:,:,ball_ind),SU_ball_alphas(:,direct_dt_ind,ball_ind));
sim_rel_alphas = bsxfun(@rdivide,sim_alphas,sim_alphas(:,sim_dt_ind));

f1 = figure();hold on
% plot_errorbar_quantiles(poss_bin_dts,Mod_rel_alphas,[25 50 75],'b');
plot_errorbar_quantiles(poss_bin_dts(use_np_bins),ball_rel_alphas(:,use_np_bins),[25 50 75],'r');
plot_errorbar_quantiles(0.01*sim_params.poss_ubins,sim_rel_alphas,[25 50 75],'k');
% legend('Nonparametric','Model-sim','Location','Northwest');
xlabel('Time window (s)');
ylabel('Relative alpha');
set(gca,'xscale','log');
xlim([0.004 4]);

%look at robustness eps-based calcs vs time window
eps_npts = arrayfun(@(x) x.eps_ball_npts(ball_ind),all_cell_data(SU_uset,:)); %number of data points at each eps value
eps_bmean = arrayfun(@(x) x.eps_ball_boot(ball_ind,1),all_cell_data(SU_uset,:)); %boot-strap mean var est
eps_bstd = arrayfun(@(x) x.eps_ball_boot(ball_ind,2),all_cell_data(SU_uset,:)); %bootstrap sd of var est
eps_cv = eps_bstd./eps_bmean; %CV of bootstrap resampling of var ests

% f2 = figure(); hold on
% % plot_errorbar_quantiles(poss_bin_dts,eps_npts,[25 50 75]);
% plot_errorbar_quantiles(poss_bin_dts,eps_cv,[25 50 75],'r');
% xlabel('Time window (s)');
% ylabel('CV');
% set(gca,'xscale','log');
% xlim([0.004 2.5]);

% 
fig_width = 5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'alpha_vs_timewin_compare.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'epsCV_vs_timewin.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% compare rate variance captured by the model with direct estimates
close all
%fraction of firing rate variance of the model compared to nonpar est
mod_R2 = Mod_tot_vars(:,mod_dt_ind)./SU_ball_vars(:,direct_dt_ind,ball_ind); 

%plot distribution of model R2 values
f1 = figure();
nbins = 20; 
bin_width = range(mod_R2)/nbins;
bin_edges = linspace(0,max(mod_R2)+bin_width/2,nbins + 1);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
n = histc(mod_R2,bin_edges);
h = bar(gca,bin_cents,n(1:end-1));
set(h,'barwidth',0.8,'faceColor','k');
xlabel('Model R2');
ylabel('Number of cells');
xlim(minmax(bin_edges));
yl = ylim();
line(median(mod_R2)+[0 0],yl,'color','b');

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Mod_R2_dist.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% DIRECT ESTIMATES OF ALPHA VS RF PROPERTIES
close all
ex_unit_ids = [20 23]; %set of units for examples
ex_set = find(ismember(SU_uset,ex_unit_ids));

%plot alpha vs RF ecc
f1 = figure()
subplot(2,2,1); hold on
plot(RF_ecc_avg,SU_ball_alphas(:,direct_dt_ind,ball_ind),'.','markersize',mSize)
[a,b] = corr(RF_ecc_avg,SU_ball_alphas(:,direct_dt_ind,ball_ind),'type','spearman');
if ~isempty(ex_set)
   plot(RF_ecc_avg(ex_set(1)),SU_ball_alphas(ex_set(1),direct_dt_ind,ball_ind),'ro');
   plot(RF_ecc_avg(ex_set(2)),SU_ball_alphas(ex_set(2),direct_dt_ind,ball_ind),'ko');
end
xlim([0 10.5]); ylim([0 1]);
title(sprintf('ECC corr; %.3f, p %.2g',a,b));
xlabel('Eccentricity (deg)');
ylabel('Alpha');
r = robustfit(RF_ecc_avg,SU_ball_alphas(:,direct_dt_ind,ball_ind)); %add robust linear reg line
xr = xlim(); xx = linspace(xr(1),xr(2),100);
plot(xx,r(1) + r(2)*xx,'r');

%plot alpha vs RF width (log scale)
subplot(2,2,2); hold on
plot(RF_avg_width,SU_ball_alphas(:,direct_dt_ind,ball_ind),'.','markersize',mSize)
xlim([0.045 1.2]); ylim([0 1]);
[a,b] = corr(log10(RF_avg_width),SU_ball_alphas(:,direct_dt_ind,ball_ind),'type','spearman');
if ~isempty(ex_set)
   plot(RF_avg_width(ex_set(1)),SU_ball_alphas(ex_set(1),direct_dt_ind,ball_ind),'ro');
   plot(RF_avg_width(ex_set(2)),SU_ball_alphas(ex_set(2),direct_dt_ind,ball_ind),'ko');
end
title(sprintf('Width corr; %.3f, p %.2g',a,b));
xlabel('RF width (deg)'); 
ylabel('Alpha');
r = robustfit(log10(RF_avg_width),SU_ball_alphas(:,direct_dt_ind,ball_ind));
xr = xlim(); xx = linspace(log10(xr(1)),log10(xr(2)),100);
plot(10.^xx,(r(1) + r(2)*xx),'r');
set(gca,'xscale','log');
set(gca,'xticklabel',10.^str2num(get(gca,'xticklabel'))); 

%plot RF ecc vs RF width (logscale)
subplot(2,2,3); hold on
plot(RF_ecc_avg,RF_avg_width,'.','markersize',mSize);
xlabel('RF ecc (deg)');
ylabel('RF width (deg)');
ylim([0.045 1.2]);
xlim([0 10.5]);
[a,b] = corr(RF_ecc_avg,RF_avg_width,'type','spearman');
if ~isempty(ex_set)
   plot(RF_ecc_avg(ex_set(1)),RF_avg_width(ex_set(1)),'ro');
   plot(RF_ecc_avg(ex_set(2)),RF_avg_width(ex_set(2)),'ko');
end
title(sprintf('Width-ecc corr; %.3f, p %.2g',a,b));
r = robustfit(RF_ecc_avg,log10(RF_avg_width));
xr = xlim(); xx = linspace(xr(1),xr(2),100);
plot(xx,10.^(r(1) + r(2)*xx),'r');
set(gca,'yscale','log');
set(gca,'yticklabel',10.^str2num(get(gca,'yticklabel'))); 

%plot RF phase sensitivity vs alpha
subplot(2,2,4); hold on
plot(RF_PRI,SU_ball_alphas(:,direct_dt_ind,ball_ind),'.','markersize',mSize)
[a,b] = corr(RF_PRI,SU_ball_alphas(:,direct_dt_ind,ball_ind),'type','spearman');
title(sprintf('PRI corr; %.3f, p %.2g',a,b));
xlabel('PRI');
ylabel('Alpha');
xlim([0 2.1]); ylim([0 1]);
if ~isempty(ex_set)
   plot(RF_PRI(ex_set(1)),SU_ball_alphas(ex_set(1),direct_dt_ind,ball_ind),'ro');
   plot(RF_PRI(ex_set(2)),SU_ball_alphas(ex_set(2),direct_dt_ind,ball_ind),'ko');
end
r = robustfit(RF_PRI,SU_ball_alphas(:,direct_dt_ind,ball_ind));
xr = xlim(); xx = linspace(xr(1),xr(2),100);
plot(xx,r(1) + r(2)*xx,'r');


% fig_width = 8; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Alpha_vs_RF.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
%% COMPARE DIRECT and EP-corrected FF ESTIMATES
close all

psth_FFs = arrayfun(@(x) x.psth_FF,all_cell_data(SU_uset,direct_used_dts));
ball_FFs = arrayfun(@(x) x.ball_FF(ball_ind),all_cell_data(SU_uset,direct_used_dts));

f1 = figure();
plot(psth_FFs(:,direct_dt_ind),ball_FFs(:,direct_dt_ind),'.','markersize',mSize)
line([0 2],[0 2],'color','k');
line([0 2],[1 1],'color','k','linestyle','--');
line([1 1],[0 2],'color','k','linestyle','--');
xlabel('PSTH-based FF');
ylabel('EP-corrected FF');

% cur_dt_ind = find(direct_bin_dts == 0.08);
% f2 = figure();
% plot(psth_FFs(:,cur_dt_ind),ball_FFs(:,cur_dt_ind),'.','markersize',mSize)
% line([0 2],[0 2],'color','k');
% line([0 2],[1 1],'color','k','linestyle','--');
% line([1 1],[0 2],'color','k','linestyle','--');
% xlabel('PSTH-based FF');
% ylabel('EP-corrected FF');

% 
% fig_width = 4; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Direct_FF_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% FF bias as a function of time window, model-nonpar compare
psth_FFs = arrayfun(@(x) x.psth_FF,all_cell_data(SU_uset,direct_used_dts));
ball_FFs = arrayfun(@(x) x.ball_FF(ball_ind),all_cell_data(SU_uset,direct_used_dts));
FF_bias = psth_FFs - ball_FFs; %nonparmetric estimate of FF bias

%get model-based estimate of FF bias
Mod_tot_vars = arrayfun(@(x) mean(x.mod_tot_vars_noEM),all_cell_data(SU_uset,:)); %PSTH variance of model-rates (avg over recs)
Mod_FF_bias = Mod_tot_vars.*Mod_alphas./SU_mean_cnts;

%get simulation-based FF calcs
[sim_FF,sim_rate_FF] = deal(nan(length(SU_uset),length(sim_params.poss_ubins)));
for ii = 1:length(SU_uset)
   cur_FF = all_cell_data(SU_uset(ii),sim_dt_ind).sim_data.FF_bias; %get model FF bias across time windows and EP SDs
   sim_FF(ii,:) = interp1(sim_params.poss_SDs(1:end-1),cur_FF(1:end-1,:),EP_SDs(ii)); %interpolate onto the observed EP SD for this rec (during rpt trials)
   cur_rate_FF = all_cell_data(SU_uset(ii),sim_dt_ind).sim_data.tot_vars./all_cell_data(SU_uset(ii),sim_dt_ind).sim_data.mean_rates; %get the firing rate FF (variance/mean ratio)
   sim_rate_FF(ii,:) = interp1(sim_params.poss_SDs(1:end-1),cur_rate_FF(1:end-1,:),EP_SDs(ii)); %interpolate onto observed EP SD
end

%normalize by values at base time-resolution
sim_FF_rel = bsxfun(@rdivide,sim_FF,sim_FF(:,1)); %normalize the simulated FFs by their value at base dt (10ms)
FF_rel = bsxfun(@rdivide,FF_bias,FF_bias(:,direct_dt_ind)); %normalize the direct FF bias by their value at base dt
mod_FF_rel = bsxfun(@rdivide,Mod_FF_bias,Mod_FF_bias(:,direct_dt_ind)); %normalize the direct FF bias by their value at base dt

% %plot absolute FF bias across time windows
% f1 = figure();
% plot_errorbar_quantiles(poss_bin_dts,FF_bias,[25 50 75]);
% hold on
% % plot_errorbar_quantiles(0.01*sim_params.poss_ubins,sim_FF,[25 50 75],'r');
% set(gca,'xscale','log');
% xlim([0.004 4]);
% xlabel('Time window (s)');
% ylabel('Relative FF bias');

%plot relative FF biases across time windows
f2 = figure();hold on
plot_errorbar_quantiles(poss_bin_dts(use_np_bins),FF_rel(:,use_np_bins),[25 50 75],'r');
plot_errorbar_quantiles(0.01*sim_params.poss_ubins,sim_FF_rel,[25 50 75],'k');
% errorbar(poss_bin_dts,nanmean(FF_rel),nanstd(FF_rel)/sqrt(length(SU_uset)),'r');
% errorbar(0.01*sim_params.poss_ubins,nanmean(sim_FF_rel),nanstd(sim_FF_rel)/sqrt(length(SU_uset)),'k');
set(gca,'xscale','log');
xlim([0.004 4]);
ylim([0.5 2.5]);
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
xlabel('Time window (s)');
ylabel('Relative FF bias');

f4 = figure(); hold on
plot_errorbar_quantiles(poss_bin_dts,psth_FFs,[25 50 75],'b');
plot_errorbar_quantiles(poss_bin_dts(use_np_bins),ball_FFs(:,use_np_bins),[25 50 75],'r');
set(gca,'xscale','log');
xlim([0.004 5]);
ylim([0.4 3]);
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
xlabel('Time window (s)');
ylabel('Fano Factor');



% %get the relative firing rate variance and its var/mean ratio as a fnx of
% %time window (for direct estimates)
% SU_rel_ball_var = bsxfun(@rdivide,squeeze(SU_ball_vars(:,:,ball_ind)),squeeze(SU_ball_vars(:,direct_dt_ind,ball_ind)));
% SU_rate_FF = bsxfun(@rdivide,squeeze(SU_ball_vars(:,:,ball_ind)),squeeze(SU_mean_cnts)); %divide rate variance by mean
% SU_rel_rate_FF = bsxfun(@rdivide,SU_rate_FF,SU_rate_FF(:,direct_dt_ind)); %normalize this 'rate-FF' by its baseline value

% %plot the dependence of 'rate-FF' on time window against the FF bias
% %(direct estimates
% f3 = figure(); hold on
% plot_errorbar_quantiles(poss_bin_dts,SU_rel_rate_FF,[25 50 75])
% plot_errorbar_quantiles(poss_bin_dts,FF_rel,[25 50 75],'r');
% set(gca,'xscale','log');
% xlim([0.004 1.5]);
% ylim([-0.5 2.5]);
% xlabel('Time window (s)');
% ylabel('Relative FF');

% f4 = figure(); hold on
% plot_errorbar_quantiles(poss_bin_dts,SU_rel_ball_var,[25 50 75])
% plot(poss_bin_dts,poss_bin_dts./dt_ind,'k') 
% set(gca,'xscale','log');
% xlim([0.0025 1.5]);
% 
fig_width = 5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'FF_timewin.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

% figufy(f2);
% fname = [fig_dir 'FFbias_timewin.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

% figufy(f3);
% fname = [fig_dir 'rateFF_FFbias_timewin.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

% figufy(f4);
% fname = [fig_dir 'FF_timewin.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% look at size of FF bias and alpha as a function of time window, using full model simulations
close all

%actual EP SDs across recordings
sim_EP_SDs = arrayfun(@(x) x.sim_params.poss_SDs(end),all_cell_data(SU_uset,1));
target_SD = nanmedian(sim_EP_SDs);
sim_dt_ind = find(sim_params.poss_ubins == 1);

[alpha_timefun,FF_bias_timefun] = deal(nan(length(SU_uset),length(sim_params.poss_ubins)));
for ii = 1:length(SU_uset)
   cur_FF_bias = all_cell_data(SU_uset(ii),sim_dt_ind).sim_data.FF_bias; 
   cur_alpha = all_cell_data(SU_uset(ii),sim_dt_ind).sim_data.alphas; 
   FF_bias_timefun(ii,:) = interp1(sim_params.poss_SDs(1:end-1),cur_FF_bias(1:end-1,:),target_SD);
   alpha_timefun(ii,:) = interp1(sim_params.poss_SDs(1:end-1),cur_alpha(1:end-1,:),target_SD);
end

sim_FF_rel = bsxfun(@rdivide,FF_bias_timefun,FF_bias_timefun(:,1));
sim_alpha_rel = bsxfun(@rdivide,alpha_timefun,alpha_timefun(:,1));
% sim_FF_bias = cell2mat(arrayfun(@(x) x.sim_data.FF_bias(use_SD_ind,:),all_cell_data(SU_uset,direct_dt_ind),'uniformoutput',0));
% sim_alpha = cell2mat(arrayfun(@(x) x.sim_data.alphas(use_SD_ind,:),all_cell_data(SU_uset,direct_dt_ind),'uniformoutput',0));
% sim_tot_vars = cell2mat(arrayfun(@(x) x.sim_data.tot_vars(use_SD_ind,:),all_cell_data(SU_uset,direct_dt_ind),'uniformoutput',0));

% sim_FF_rel = bsxfun(@rdivide,sim_FF_bias,sim_FF_bias(:,1));
% sim_alpha_rel = bsxfun(@rdivide,sim_alpha,sim_alpha(:,1));
% sim_totvars_rel = bsxfun(@rdivide,sim_tot_vars,sim_tot_vars(:,1));

f1 = figure();
subplot(2,1,1);
errorbar(1e3*sim_params.poss_ubins*poss_bin_dts(direct_dt_ind),nanmean(sim_FF_rel),nanstd(sim_FF_rel)/sqrt(length(SU_uset)));
set(gca,'xscale','log');
xlim([7.5 2500]);
xlabel('Time window (ms)');
ylabel('Fano Factor bias');
ylim([.75 2])

subplot(2,1,2);
errorbar(1e3*sim_params.poss_ubins*poss_bin_dts(direct_dt_ind),nanmean(sim_alpha_rel),nanstd(sim_alpha_rel)/sqrt(length(SU_uset)));
set(gca,'xscale','log');
xlim([7.5 2500]);
ylim([0.8 1.2]);
xlabel('Time window (ms)');
ylabel('Relative Alpha');

% f3 = figure();
% errorbar(1e3*sim_params.poss_ubins*poss_bin_dts(direct_dt_ind),nanmean(sim_totvars_rel),nanstd(sim_totvars_rel));
% hold on
% plot(1e3*sim_params.poss_ubins*poss_bin_dts(direct_dt_ind),sim_params.poss_ubins.^2,'r')
% set(gca,'xscale','log');
% xlim([5 1500]);
% xlabel('Time window (ms)');
% ylabel('Relative rate variance');

% fig_width = 4; rel_height = 1.6;
% figufy(f1);
% fname = [fig_dir 'Modsim_timebin_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% look at FF bias and alpha as a function of EP SD
close all

sim_FF_bias = cell2mat(arrayfun(@(x) x.sim_data.FF_bias(:,sim_dt_ind)',all_cell_data(SU_uset,1),'uniformoutput',0));
sim_alpha = cell2mat(arrayfun(@(x) x.sim_data.alphas(:,sim_dt_ind)',all_cell_data(SU_uset,1),'uniformoutput',0));

%normalize by the median EP SD 
norm_SD = nanmedian(EP_SDs);
sim_alpha_norm = interp1(sim_params.poss_SDs(1:end-1),sim_alpha(:,1:end-1)',norm_SD);
sim_FF_norm = interp1(sim_params.poss_SDs(1:end-1),sim_FF_bias(:,1:end-1)',norm_SD);
sim_FF_rel = bsxfun(@rdivide,sim_FF_bias(:,1:end-1),sim_FF_norm');
sim_alpha_rel = bsxfun(@rdivide,sim_alpha(:,1:end-1),sim_alpha_norm');


observed_EP_range = minmax(EP_SDs);
interp_alphas = interp1(sim_params.poss_SDs(1:end-1),sim_alpha(:,1:end-1)',observed_EP_range');
% EPrange_fold_change = diff(interp_alphas)./interp_alphas(1,:);
EPrange_fold_change = bsxfun(@rdivide,interp_alphas,sim_alpha_norm);

sd_range = [0.0 0.22]; %range for plotting

%plot the relative change in EM-induced tbt variance with EP SD
f1 = figure(); 
% errorbar(sim_params.poss_SDs(1:end-1),nanmean(sim_alpha_rel),nanstd(sim_alpha_rel),'k');
plot_errorbar_quantiles(sim_params.poss_SDs(1:end-1),sim_alpha_rel,[25 50 75],'k');
xlabel('Eye position SD (deg)');
ylabel('Relative EM-induced variance');
xlim(sd_range);
ylim([0. 1.6])

%plot the distribution of EP SDs
f2 = figure();
nbins = 15;
bin_width = range(EP_SDs)/nbins;
bin_edges = linspace(min(EP_SDs)-bin_width/2,max(EP_SDs)+bin_width/2,nbins + 1);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
n = histc(EP_SDs,bin_edges);
h = bar(gca,bin_cents,n(1:end-1));
set(h,'barwidth',0.8,'faceColor','k');
xlabel('Eye position SD (deg)');
xlim(minmax(bin_edges));
yl = ylim();
line(median(EP_SDs)+[0 0],yl,'color','b');
xlim(sd_range);

% 
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Modsim_EPSD_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 4; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'EPSD_dist.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);


%% noise correlation analysis
% close all
dt_ind = 0.01; %which dt value to use
direct_dt_ind = find(direct_bin_dts == dt_ind);

norm_type = 'fixed'; %normalize by same value (direct, psth-based noise var estimates)
% norm_type = 'noise'; %normalize by individual estimates of noise variances

dt = direct_bin_dts(direct_dt_ind);
tlags = all_cell_data(SU_uset(1),direct_dt_ind).tlags;
cent_lag = find(tlags == 0);

%get properties of each SU pair
SU_CID = [all_cell_data(SU_uset).cell_ID];
[pair_RF_eccs,pair_RF_widths,pair_SU_chs,pair_psth_Nvars,pair_ball_Nvars] = deal(nan(length(upairs),2));
for ii = 1:length(upairs)
    curset = find(ismember(SU_CID,all_pair_data(upairs(ii),1).cell_IDs));
    pair_RF_eccs(ii,:) = RF_ecc_avg(curset);
    pair_RF_widths(ii,:) = RF_avg_width(curset);
    pair_SU_chs(ii,1) = all_cell_data(SU_uset(curset(1)),1).unit_data.probe_number;
    pair_SU_chs(ii,2) = all_cell_data(SU_uset(curset(2)),1).unit_data.probe_number;    
    pair_psth_Nvars(ii,:) = SU_psth_noisevars(curset,direct_dt_ind); %psth-based noise variance estimates
    pair_ball_Nvars(ii,:) = SU_ball_noisevars(curset,direct_dt_ind,ball_ind); %EP-corrected noise variance estimates
end
pair_rpt_trials = arrayfun(@(x) length(x.pair_rpt_set),all_pair_data(upairs,1)); %number of rpt trials with both units isolated
% pair_exptnum = arrayfun(@(x) x.Expt_num,all_pair_data(upairs,1));

%selection criteria for pairs, 
cur_upair_inds = find(pair_rpt_trials > min_nTrials); %make sure both units had at least Ntrials simultaneous rec
cur_upairs = upairs(cur_upair_inds); 
fprintf('Using %d pairs\n',length(cur_upairs));

%get mean RF widths
arith_mean_width = mean(pair_RF_widths,2);
% arith_mean_width = mean(log10(pair_RF_widths),2);
arith_mean_width = arith_mean_width(cur_upair_inds);

%total pairwise cross covariance function
all_tot_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.tot_xcovar,1),all_pair_data(cur_upairs,direct_dt_ind),'uniformoutput',0)));
%get psth-based and EP-based cross-covars
all_psth_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.pair_xcovar,1),all_pair_data(cur_upairs,direct_dt_ind),'uniformoutput',0)));
all_EP_xcovs = squeeze(cat(1,cell2mat(arrayfun(@(x) mean(x.eps_xcovar_LOO(:,ball_ind,:),1),all_pair_data(cur_upairs,direct_dt_ind),'uniformoutput',0))));

%get variance normalization (product of PSTH-based noise variance
%estimates)
all_xcov_norms = cat(1,cell2mat(arrayfun(@(x) mean(x.at_var_norm),all_pair_data(cur_upairs,direct_dt_ind),'uniformoutput',0)));
all_xcov_pnorms = sqrt(prod(pair_psth_Nvars(cur_upair_inds,:),2)); %noise-variance normalization based on PSTH estimate
all_xcov_enorms = sqrt(prod(pair_ball_Nvars(cur_upair_inds,:),2)); %noise-variance normalization based on EP-corrected estimate

%normalize to correlations
if strcmp(norm_type,'fixed') %normalize by same noise-var estimate
    all_EP_xcorrs = bsxfun(@rdivide,all_EP_xcovs,all_xcov_norms);
    all_psth_xcorrs = bsxfun(@rdivide,all_psth_xcovs,all_xcov_norms);
    %compute noise correlations (PSTH and EP-based)
    all_EP_noisecorrs = bsxfun(@rdivide,(all_tot_xcovs - all_EP_xcovs),all_xcov_norms);
    all_psth_noisecorrs = bsxfun(@rdivide,(all_tot_xcovs - all_psth_xcovs),all_xcov_norms);
elseif strcmp(norm_type,'noise') %normalize by individual noise-variance estimates
    all_EP_xcorrs = bsxfun(@rdivide,all_EP_xcovs,all_xcov_enorms);
    all_psth_xcorrs = bsxfun(@rdivide,all_psth_xcovs,all_xcov_pnorms);
    %compute noise correlations (PSTH and EP-based)
    all_EP_noisecorrs = bsxfun(@rdivide,(all_tot_xcovs - all_EP_xcovs),all_xcov_enorms);
    all_psth_noisecorrs = bsxfun(@rdivide,(all_tot_xcovs - all_psth_xcovs),all_xcov_pnorms);
else
    error('invalid norm type');
end

%grab correlation-values at central time lag
EP_xcorrs_cent = squeeze(all_EP_xcorrs(:,cent_lag));
psth_xcorrs_cent = squeeze(all_psth_xcorrs(:,cent_lag));
EP_noisecorrs_cent = squeeze(all_EP_noisecorrs(:,cent_lag));
psth_noisecorrs_cent = squeeze(all_psth_noisecorrs(:,cent_lag));

%get linear fits to extract fraction of covariances captured by PSTH
[psth_fract,EP_fract,psth_EP_fract] = deal(nan(length(cur_upairs),1));
for ii = 1:length(cur_upairs)
    if length(tlags) > 1
        B = regress(all_psth_xcovs(ii,:)',[all_tot_xcovs(ii,:)' ones(length(tlags),1)]);
        psth_fract(ii) = B(1); %fraction of total covariance captured by PSTH
        B = regress(all_EP_xcovs(ii,:)',[ all_tot_xcovs(ii,:)' ones(length(tlags),1)]);
        EP_fract(ii) = B(1); %fraction of total covariance captured by EP-corrected
        B = regress(all_psth_xcovs(ii,:)',[all_EP_xcovs(ii,:)' ones(length(tlags),1)]);
        psth_EP_fract(ii) = B(1); %fraction of EP-corrected captured by PSTH
    end
end

%get signal-noise correlations
[EP_sig_noise_corr,EP_p] = corr(EP_noisecorrs_cent,EP_xcorrs_cent,'type','spearman');
[PSTH_sig_noise_corr,PSTH_p] = corr(psth_noisecorrs_cent,psth_xcorrs_cent,'type','spearman');
fprintf('PSTH cent corr: %.3f p: %.3f\n',PSTH_sig_noise_corr,PSTH_p);
fprintf('EP cent corr: %.3f p: %.3f\n',EP_sig_noise_corr,EP_p);

xl1 = [-0.3 0.5]; 
yl1 = [-0.3 0.3];
xx = linspace(-0.3,0.5,100);
nboots = 5000;
%plot sig/noise corr relationships
f1 = figure(); 
subplot(2,1,1) %first for PSTH-based
hold on
boot_ypred = nan(nboots,length(xx));
n_dpts = length(psth_xcorrs_cent);
for ii = 1:nboots 
    boot_samp = randi(n_dpts,n_dpts,1);
    r1 = regress(psth_noisecorrs_cent(boot_samp),[ones(n_dpts,1) psth_xcorrs_cent(boot_samp)]);
    boot_ypred(ii,:) = r1(1) + r1(2)*xx;
end
% r1 = regress(psth_noisecorrs_cent,[ones(size(psth_xcorrs_cent)) psth_xcorrs_cent]);
% plot(xx,r1(1) + r1(2)*xx,'r--')
UE = prctile(boot_ypred,97.5) - mean(boot_ypred);
LE = mean(boot_ypred) - prctile(boot_ypred,2.5);
shadedErrorBar(xx,mean(boot_ypred),[UE; LE]);
plot(psth_xcorrs_cent,psth_noisecorrs_cent,'r.','markersize',mSize);
line(xl1,[0 0],'color','k','linestyle','--'); line([0 0],yl1,'color','k','linestyle','--');
% line([-0.5 0.5],[-0.5 0.5],'color','k','linestyle','--');
% r1 = robustfit(psth_xcorrs_cent,psth_noisecorrs_cent);
xlim(xl1); ylim(yl1);
xlabel('Signal correlation');
ylabel('Noise correlation');

subplot(2,1,2) %now for EP-corrected
hold on
boot_ypred = nan(nboots,length(xx));
n_dpts = length(EP_noisecorrs_cent);
for ii = 1:nboots
    boot_samp = randi(n_dpts,n_dpts,1);
    r1 = regress(EP_noisecorrs_cent(boot_samp),[ones(n_dpts,1) EP_xcorrs_cent(boot_samp)]);
    boot_ypred(ii,:) = r1(1) + r1(2)*xx;
end
% r1 = regress(psth_noisecorrs_cent,[ones(size(psth_xcorrs_cent)) psth_xcorrs_cent]);
% plot(xx,r1(1) + r1(2)*xx,'r--')
UE = prctile(boot_ypred,97.5) - mean(boot_ypred);
LE = mean(boot_ypred) - prctile(boot_ypred,2.5);
shadedErrorBar(xx,mean(boot_ypred),[UE; LE]);
plot(EP_xcorrs_cent,EP_noisecorrs_cent,'b.','markersize',mSize);
line(xl1,[0 0],'color','k','linestyle','--'); line([0 0],yl1,'color','k','linestyle','--');
% line([-0.5 0.5],[-0.5 0.5],'color','k','linestyle','--');
% r1 = robustfit(EP_xcorrs_cent,EP_noisecorrs_cent);
% plot(xx,r1(1) + r1(2)*xx,'b--')
xlim(xl1); ylim(yl1);
xlabel('Signal correlation');
ylabel('Noise correlation');

% %plot EP vs psth sig corr estimates
% f2 = figure();
% plot(EP_xcorrs_cent,psth_xcorrs_cent,'.','markersize',mSize);
% r = robustfit(EP_xcorrs_cent,psth_xcorrs_cent);
% hold on
% plot(xx,r(1) + r(2)*xx,'r');
% xlim(xl1); ylim(xl1);
% line(xl1,xl1,'color','k');
% xlabel('EP-corrected signal correlation');
% ylabel('PSTH-based signal correlation');
% 
%plot 'alpha' vs mean RF width
if length(tlags) > 1
    f3 = figure();
    plot(arith_mean_width,psth_EP_fract,'.','markersize',mSize);
    hold on
%     [~,ord] = sort(arith_mean_width);
%     smooth_bnd = 80;
%     plot(arith_mean_width(ord),smooth(psth_EP_fract(ord),smooth_bnd,'lowess'),'r');
    xlim([0.05 1])
    r = robustfit(log10(arith_mean_width),psth_EP_fract);
    xr = xlim();
    xx = linspace(log10(xr(1)),log10(xr(2)),100);
    plot(10.^xx,r(1) + r(2)*xx,'r');
    xlabel('Mean RF width (deg)');
    ylabel('Fraction signal variance');
    set(gca,'xscale','log');
    set(gca,'xtickLabel',10.^str2num(get(gca,'xticklabel')));
end

%plot the fraction of the total xcov function captured by the psth- and
%EP-based estimators of rate covariance
if length(tlags) > 1
    f4 = figure();
    plot(psth_fract,EP_fract,'.','markersize',mSize);
    xl = xlim(); yl = ylim();
    line(yl,yl,'color','k');
    xlabel('PSTH variance fraction');
    ylabel('EP-corrected variance fraction');
end
% 
% fig_width = 4; rel_height = 2;
% figufy(f1);
% fname = [fig_dir 'signoise_Xcorr.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% % fig_width = 4; rel_height = 1;
% % figufy(f2);
% % fname = [fig_dir 'psth_EP_sigcorr_compare.pdf'];
% % exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f2);
% % 
% fig_width = 4; rel_height = 1;
% figufy(f3);
% fname = [fig_dir 'xcorr_alpha_RFwidth.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

% fig_width = 4; rel_height = 1;
% figufy(f4);
% fname = [fig_dir 'EP_PSTH_varfracts.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);
% 
%% plot noise cov frac vs time window
close all

% norm_type = 'fixed'; %normalize by same value (direct, psth-based noise var estimates)
norm_type = 'noise'; %normalize by individual estimates of noise variances
n_boot_samps = 50; %number of boostrap samples for estimating error bars in rho

%get SU pair stats
SU_CID = [all_cell_data(SU_uset).cell_ID];
[pair_RF_eccs,pair_RF_widths,pair_SU_chs,pair_EPSDs] = deal(nan(length(upairs),2));
[pair_psth_Nvars,pair_ball_Nvars] = deal(nan(length(upairs),2,length(direct_bin_dts)));
for ii = 1:length(upairs)
    curset = find(ismember(SU_CID,all_pair_data(upairs(ii),1).cell_IDs));
    pair_RF_eccs(ii,:) = RF_ecc_avg(curset);
    pair_RF_widths(ii,:) = RF_avg_width(curset);
    pair_EPSDs(ii,:) = EP_SDs(curset);
    pair_SU_chs(ii,1) = all_cell_data(SU_uset(curset(1)),1).unit_data.probe_number;
    pair_SU_chs(ii,2) = all_cell_data(SU_uset(curset(2)),1).unit_data.probe_number;
    pair_psth_Nvars(ii,:,:) = SU_psth_noisevars(curset,:); %psth-based noise variance estimate
    pair_ball_Nvars(ii,:,:) = SU_ball_noisevars(curset,:,ball_ind); %EP-corrected noise variance estimates
end
pair_rpt_trials = arrayfun(@(x) length(x.pair_rpt_set),all_pair_data(upairs,1));
pair_exptnum = arrayfun(@(x) x.Expt_num,all_pair_data(upairs,1));

%selection criteria for pairs, 
cur_upair_inds = find(pair_rpt_trials > min_nTrials); %make sure both units had at least Ntrials simultaneous rec
cur_upairs = upairs(cur_upair_inds); 

psth_paircov_frac = nan(length(direct_bin_dts),1);
n_bad_norms = nan(length(direct_bin_dts),1);
[psth_sig_noise_rho,EP_sig_noise_rho] = deal(nan(length(direct_bin_dts),6));
for tt = 1:length(direct_bin_dts)
    
    tlags = all_cell_data(SU_uset(1),tt).tlags; %time lag axis at this time-res
    cent_lag = find(tlags == 0); %central time bin
    
    %total pairwise cross covariance functio
    all_tot_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.tot_xcovar,1),all_pair_data(cur_upairs,tt),'uniformoutput',0)));
    %get psth-based and EP-based cross-covars
    all_psth_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.pair_xcovar,1),all_pair_data(cur_upairs,tt),'uniformoutput',0)));
    all_EP_xcovs = squeeze(cat(1,cell2mat(arrayfun(@(x) mean(x.eps_xcovar_LOO(:,ball_ind,:),1),all_pair_data(cur_upairs,tt),'uniformoutput',0))));
    
    %noise covs by subtracting rate-cov estimates from total xcov
    all_psth_noisecovs = all_tot_xcovs - all_psth_xcovs;
    all_EP_noisecovs = all_tot_xcovs - all_EP_xcovs;
    
    %get variance normalization (product of PSTH-based noise variance
    %estimates)
    all_xcov_norms = cat(1,cell2mat(arrayfun(@(x) mean(x.at_var_norm),all_pair_data(cur_upairs,tt),'uniformoutput',0)));
    all_xcov_pnorms = squeeze(sqrt(prod(pair_psth_Nvars(cur_upair_inds,:,tt),2)));
    all_xcov_enorms = squeeze(sqrt(prod(pair_ball_Nvars(cur_upair_inds,:,tt),2)));
    all_xcov_enorms(imag(all_xcov_enorms) > 0) = nan;
    if strcmp(norm_type,'fixed')
        all_psth_noisecorrs = bsxfun(@rdivide,all_psth_noisecovs,all_xcov_norms);
        all_EP_noisecorrs = bsxfun(@rdivide,all_EP_noisecovs,all_xcov_norms);
        all_psth_sigcorrs = bsxfun(@rdivide,all_psth_xcovs,all_xcov_norms);
        all_EP_sigcorrs = bsxfun(@rdivide,all_EP_xcovs,all_xcov_norms);
    elseif strcmp(norm_type,'noise')
        all_psth_noisecorrs = bsxfun(@rdivide,all_psth_noisecovs,all_xcov_pnorms);
        all_EP_noisecorrs = bsxfun(@rdivide,all_EP_noisecovs,all_xcov_enorms);
        all_psth_sigcorrs = bsxfun(@rdivide,all_psth_xcovs,all_xcov_pnorms);
        all_EP_sigcorrs = bsxfun(@rdivide,all_EP_xcovs,all_xcov_enorms);
    else
        error('invalid norm type');
    end
    
    avg_psth_noisecorrs(tt) = mean(all_psth_noisecorrs(:,cent_lag));
    avg_EP_noisecorrs(tt) = mean(all_EP_noisecorrs(:,cent_lag));
    med_psth_noisecorrs(tt) = median(all_psth_noisecorrs(:,cent_lag));
    med_EP_noisecorrs(tt) = median(all_EP_noisecorrs(:,cent_lag));
    
    if sum(~isnan(all_EP_xcovs(:,cent_lag))) > 0
        %regression slope betwen PSTH and EP-based xcov estimates
        r = robustfit(all_EP_xcovs(:,cent_lag),all_psth_xcovs(:,cent_lag));
        psth_paircov_frac(tt) = r(2);
        
        %get boostramp estimates of the relationship between signal and noise
        %correlation
        bootstat = bootstrp(n_boot_samps,@(x) corr(x(:,1),x(:,2),'type','spearman'),[all_psth_noisecorrs(:,cent_lag),all_psth_sigcorrs(:,cent_lag)]);
        psth_sig_noise_rho(tt,1) = nanmean(bootstat);
        psth_sig_noise_rho(tt,2) = nanstd(bootstat);
        [psth_sig_noise_rho(tt,3),psth_sig_noise_rho(tt,4)] = corr(all_psth_noisecorrs(:,cent_lag),all_psth_sigcorrs(:,cent_lag),'type','spearman');
        psth_sig_noise_rho(tt,5:6) = prctile(bootstat,[2.5 97.5]);
        
        bootstat = bootstrp(n_boot_samps,@(x) corr(x(:,1),x(:,2),'type','spearman'),[all_EP_noisecorrs(:,cent_lag),all_EP_sigcorrs(:,cent_lag)]);
        EP_sig_noise_rho(tt,1) = nanmean(bootstat);
        EP_sig_noise_rho(tt,2) = nanstd(bootstat);
        [EP_sig_noise_rho(tt,3),EP_sig_noise_rho(tt,4)] = corr(all_EP_noisecorrs(:,cent_lag),all_EP_sigcorrs(:,cent_lag),'type','spearman');
        EP_sig_noise_rho(tt,5:6) = prctile(bootstat,[2.5 97.5]);
        
        n_bad_norms(tt) = sum(isnan(all_xcov_enorms));
    end
end

%get signal-noise corr relationship for simulations
sim_psth_sig_noise_rho = nan(length(sim_params.poss_ubins),6);
for tt = 1:length(sim_params.poss_ubins)
    all_tot_xcovs = cell2mat(arrayfun(@(x) x.sim_data.tot_xcovar(:,tt)',all_pair_data(cur_upairs,1),'uniformoutput',0));
    all_psth_xcovs = cell2mat(arrayfun(@(x) x.sim_data.psth_xcovar(:,tt)',all_pair_data(cur_upairs,1),'uniformoutput',0));
    all_covar_norm = cell2mat(arrayfun(@(x) x.sim_data.covar_norm(:,tt)',all_pair_data(cur_upairs,1),'uniformoutput',0));
    
    %interpolate onto EP SD (avg of pair) from repeat trials
    cur_pair_EP = mean(pair_EPSDs(cur_upair_inds,:),2); %avg rpt EP SD of SU pair
    [cur_tot_xcovs,cur_psth_xcovs,cur_covar_norm] = deal(nan(length(cur_upairs),1));
    for ii = 1:length(cur_upair_inds)
        cur_tot_xcovs(ii) = interp1(sim_params.poss_SDs(1:end-1)',all_tot_xcovs(ii,1:end-1),cur_pair_EP(ii));
        cur_psth_xcovs(ii) = interp1(sim_params.poss_SDs(1:end-1)',all_psth_xcovs(ii,1:end-1),cur_pair_EP(ii));
        cur_covar_norm(ii) = interp1(sim_params.poss_SDs(1:end-1)',all_covar_norm(ii,1:end-1),cur_pair_EP(ii));
    end
    cur_psth_noisecovs = cur_tot_xcovs - cur_psth_xcovs; %psth-based noise cov estimates
    
    %normalize by across-trial variance
    sim_psth_noisecorrs = cur_psth_noisecovs./cur_covar_norm;
    sim_psth_sigcorrs = cur_psth_xcovs./cur_covar_norm;
    
    %get boostramp estimates of the relationship between signal and noise
    %correlation
    bootstat = bootstrp(n_boot_samps,@(x) corr(x(:,1),x(:,2),'type','spearman'),[sim_psth_noisecorrs,sim_psth_sigcorrs]);
    sim_psth_sig_noise_rho(tt,1) = nanmean(bootstat);
    sim_psth_sig_noise_rho(tt,2) = nanstd(bootstat);
    [sim_psth_sig_noise_rho(tt,3),sim_psth_sig_noise_rho(tt,4)] = corr(sim_psth_noisecorrs,sim_psth_sigcorrs,'type','spearman');
    sim_psth_sig_noise_rho(tt,5:6) = prctile(bootstat,[2.5 97.5]);
    
end

use_direct_dt_range = find(direct_bin_dts < 0.25); %set of time windows to use for this analysis

%plot sig-nois correlations as a function of time window with 95%CI bars
f1 = figure();hold on
errorbar(direct_bin_dts(use_direct_dt_range),psth_sig_noise_rho(use_direct_dt_range,1),...
    psth_sig_noise_rho(use_direct_dt_range,1) - psth_sig_noise_rho(use_direct_dt_range,5),...
    psth_sig_noise_rho(use_direct_dt_range,6) - psth_sig_noise_rho(use_direct_dt_range,1),'b');
plot(direct_bin_dts(use_direct_dt_range),psth_sig_noise_rho(use_direct_dt_range,1),'bo-');
errorbar(direct_bin_dts(use_direct_dt_range),EP_sig_noise_rho(use_direct_dt_range,1),...
    EP_sig_noise_rho(use_direct_dt_range,1) - EP_sig_noise_rho(use_direct_dt_range,5),...
    EP_sig_noise_rho(use_direct_dt_range,6) - EP_sig_noise_rho(use_direct_dt_range,1),'r')
plot(direct_bin_dts(use_direct_dt_range),EP_sig_noise_rho(use_direct_dt_range,1),'ro-');
errorbar(sim_params.poss_ubins*dt,sim_psth_sig_noise_rho(:,1),sim_psth_sig_noise_rho(:,1) - sim_psth_sig_noise_rho(:,5),...
    sim_psth_sig_noise_rho(:,6) - sim_psth_sig_noise_rho(:,1),'k')
plot(sim_params.poss_ubins*dt,sim_psth_sig_noise_rho(:,1),'ko-');
xlim([0.004 5]);
xl = xlim(); set(gca,'xscale','log');
line(xl,[0 0],'color','k','linestyle','--');
xlabel('Time window (s)');
ylabel('Signal-noise correlation');


% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'EP_PSTH_signoise_timewin.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
% sum_lags = find(abs(tlags) <= 0);
% mod_tot_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.mod_tot_covar,1),all_pair_data(upairs,mod_dt_ind),'uniformoutput',0)));
% mod_psth_xcovs = cat(1,cell2mat(arrayfun(@(x) mean(x.mod_psth_covar,1),all_pair_data(upairs,mod_dt_ind),'uniformoutput',0)));
% mod_noise_xcovs = mod_tot_xcovs - mod_psth_xcovs;
% mod_noisevar_norms = nan(length(upairs),1);
% for ii = 1:length(upairs)
%     curset = SU_uset(ismember(SU_CID,all_pair_data(upairs(ii),1).cell_IDs));
%     noise_var_1 = mean(all_cell_data(curset(1),1).mod_ep_vars) + all_cell_data(curset(1),1).ov_avg_BS; 
%     noise_var_2 = mean(all_cell_data(curset(2),1).mod_ep_vars) + all_cell_data(curset(2),1).ov_avg_BS; 
%     mod_noisevar_norms(ii) = sqrt(noise_var_1*noise_var_2);
% end
% mod_tot_xcorrs = bsxfun(@rdivide,mod_tot_xcovs,mod_noisevar_norms);
% mod_psth_xcorrs = bsxfun(@rdivide,mod_psth_xcovs,mod_noisevar_norms);

% f1 = figure(); hold on
% plot(mean(mod_tot_xcorrs(:,sum_lags),2),mean(mod_psth_xcorrs(:,sum_lags),2),'.');
% xlim(xl1); ylim(yl1);
% r = robustfit(mean(mod_tot_xcorrs(:,sum_lags),2),mean(mod_psth_xcorrs(:,sum_lags),2))
% line(xl1,yl1,'color','k');
% plot(xx,r(1) + r(2)*xx,'r')
% 
poss_ubins = sim_params.poss_ubins;
poss_SDs = sim_params.poss_SDs;
mod_tot_xcovs = cell2mat(arrayfun(@(x) reshape(x.sim_data.tot_xcovar,1,length(poss_SDs),length(poss_ubins)),all_pair_data(upairs,mod_dt_ind),'uniformoutput',0));
mod_psth_xcovs = cell2mat(arrayfun(@(x) reshape(x.sim_data.psth_xcovar,1,length(poss_SDs),length(poss_ubins)),all_pair_data(upairs,mod_dt_ind),'uniformoutput',0));
mod_noise_xcovs = mod_tot_xcovs - mod_psth_xcovs;
mod_xcov_norms = cell2mat(arrayfun(@(x) reshape(x.sim_data.covar_norm,1,length(poss_SDs),length(poss_ubins)),all_pair_data(upairs,mod_dt_ind),'uniformoutput',0));

target_SD = 0.125;
SD_ind = find(poss_SDs == target_SD);
ubin_ind = find(poss_ubins == 1);

mod_tot_xcorrs = squeeze(mod_tot_xcovs(:,SD_ind,ubin_ind)./mod_xcov_norms(:,SD_ind,ubin_ind));
mod_psth_xcorrs = squeeze(mod_psth_xcovs(:,SD_ind,ubin_ind)./mod_xcov_norms(:,SD_ind,ubin_ind));
mod_psth_noisecorrs = squeeze((mod_tot_xcovs(:,SD_ind,ubin_ind) - mod_psth_xcovs(:,SD_ind,ubin_ind))./mod_xcov_norms(:,SD_ind,ubin_ind));
mod_tot_xcovs = squeeze(mod_tot_xcovs(:,SD_ind,ubin_ind));
mod_psth_xcovs = squeeze(mod_psth_xcovs(:,SD_ind,ubin_ind));

geom_mean_width = sqrt(prod(pair_RF_widths,2));
% select_set = find(all(pair_RF_widths < 0.25,2));
select_set = 1:length(upairs);
xl1 = [-0.25 0.25]; 
yl1 = [-0.25 0.25];
xx = linspace(-0.3,0.3,100);

f1 = figure(); hold on
plot(mod_tot_xcorrs(select_set),mod_psth_xcorrs(select_set),'.');
xlim(xl1); ylim(yl1);
r = robustfit(mod_tot_xcorrs(select_set),mod_psth_xcorrs(select_set));
line(xl1,yl1,'color','k');
plot(xx,r(1) + r(2)*xx,'r')

data_use_sub = find(ismember(upairs,cur_upairs));
nbins = 6;
bin_edges = prctile(geom_mean_width,linspace(0,100,nbins+1));
[n,binids] = histc(geom_mean_width,bin_edges);
binned_r_mod = nan(nbins,1);
binned_r_data = nan(nbins,1);
avg_gm_width_mod = nan(nbins,1);
avg_gm_width_data = nan(nbins,1);
for ii = 1:nbins
    curset = find(binids == ii);
    r = robustfit(mod_tot_xcorrs(curset),mod_psth_xcorrs(curset));
    binned_r_mod(ii) = r(2);
    avg_gm_width_mod(ii) = mean(geom_mean_width(curset));

    curset = find(binids(data_use_sub) == ii);
    r = robustfit(all_EP_xcorrs_cent(curset),all_psth_xcorrs_cent(curset));
    binned_r_data(ii) = r(2);
    avg_gm_width_data(ii) = mean(geom_mean_width(data_use_sub(curset)));
end

f2 = figure();
plot(avg_gm_width_mod,binned_r_mod,'o-');
hold on
plot(avg_gm_width_data,binned_r_data,'ro-');

%%
dt = direct_bin_dts(direct_dt_ind);
close all
pair_id = 290; %[290 308 180 183]
pairloc = find(upairs == pair_id);

f1 = figure();
subplot(2,1,1); hold on
plot(tlags*dt*1e3,all_tot_xcorrs(pairloc,:),'ko-');
plot(tlags*dt*1e3,all_psth_xcorrs(pairloc,:),'bo-');
plot(tlags*dt*1e3,all_psth_noisecorrs(pairloc,:),'ro-');
ylabel('Correlation');
xlabel('Time lag (ms)');
yl = ylim();
if pair_id == 290
    yl = [-0.2 0.4];
end
ylim(yl);

subplot(2,1,2); hold on
plot(tlags*dt*1e3,all_tot_xcorrs(pairloc,:),'ko-');
plot(tlags*dt*1e3,all_EP_xcorrs(pairloc,:),'bo-');
plot(tlags*dt*1e3,all_EP_noisecorrs(pairloc,:),'ro-');
ylim(yl);
ylabel('Correlation');
xlabel('Time lag (ms)');



% fig_width = 4; rel_height = 2;
% figufy(f1);
% fname = [fig_dir sprintf('Xcorr_example_%d.pdf',pair_id)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
%%
% dt = direct_bin_dts(direct_dt_ind);
% close all
% 
% % curset = find(pair_exptnum == 10);
% curset = [17 75 135 170];
% for ii = 1:length(curset)
%     upairs(curset(ii))
%     all_pair_data(upairs(curset(ii))).Expt_num
%     subplot(2,1,1); hold on
%     plot(tlags*dt*1e3,all_psth_xcovs(curset(ii),:));
%     plot(tlags*dt*1e3,all_tot_xcovs(curset(ii),:),'k');
%     plot(tlags*dt*1e3,all_tot_xcovs(curset(ii),:)-all_psth_xcovs(curset(ii),:),'r')
%     
%     subplot(2,1,2); hold on
%     plot(tlags*dt*1e3,all_EP_xcovs(curset(ii),:));
%     plot(tlags*dt*1e3,all_tot_xcovs(curset(ii),:),'k');
%     plot(tlags*dt*1e3,all_tot_xcovs(curset(ii),:)-all_EP_xcovs(curset(ii),:),'r')
%     pause
%     clf
% end

%%
% close all
% f1 = figure();
% f2 = figure();
% for ss = 1:length(SU_uset)
%     cur_cell = all_cell_data(SU_uset(ss),direct_dt_ind);
%     figure(f1); clf;
%     plot_NMM_filters_1d(cur_cell.bestGQM,[],[],[],f1);
%     
%     figure(f2); clf; hold on
%     plot(cur_cell.EP_bin_centers,cur_cell.var_ep_binned,'.');
%     seval = cur_cell.spline_looEP.evalAt(cur_cell.eval_xx);
%     plot(cur_cell.eval_xx,seval,'r');
%     
%     pause
% end
% 
%% plot example model fits 
close all
% f1 = figure(); 
% f2 = figure();
% % f3 = figure();
% for ii = 12:length(SU_uset)
%     SU_exptNumber(ii)
%     figure(f1); clf
%     hold on
%     plot(Mod_alphas(:,mod_dt_ind),SU_ball_alphas(:,direct_dt_ind,ball_ind),'.','markersize',mSize);
%     plot(Mod_alphas(ii,mod_dt_ind),SU_ball_alphas(ii,direct_dt_ind,ball_ind),'ro','markersize',mSize);
%     line([0 1],[0 1]);
% 
%     npix = all_cell_data(SU_uset(ii),1).modFitParams.use_nPix_us;
%     pix_dx = all_cell_data(SU_uset(ii),1).modFitParams.sp_dx;
%     pix_ax = (1:npix)*pix_dx; pix_ax = pix_ax - mean(pix_ax);
%     flen = all_cell_data(SU_uset(ii),1).modFitParams.flen;
%     dt = all_cell_data(SU_uset(ii),1).modFitParams.dt;
%     tax = (1:flen)*dt - dt/2;
%     
%     figure(f2);clf
%     plot_NMM_filters_1d(all_cell_data(SU_uset(ii),1).bestGQM,pix_ax,tax,[],f2);
%     
%     [ii RF_avg_width(ii) RF_PRI(ii) RF_ecc_avg(ii)]
%     pause
% end  

ex_unit_ids = [20 23]; %set of units for examples
xr1 = [-0.5 0.5]; %range of pixel axis
xr2 = [-0.25 0.25]; %range of pixel axis
tr = [0 0.12]; %range of time axis
ex_set = find(ismember(SU_uset,ex_unit_ids));

%plot model vs direct alphas, and highlight example units
figure(f1); clf
hold on
plot(Mod_alphas(:,mod_dt_ind),SU_ball_alphas(:,direct_dt_ind,ball_ind),'.','markersize',mSize);
plot(Mod_alphas(ex_set(1),mod_dt_ind),SU_ball_alphas(ex_set(1),direct_dt_ind,ball_ind),'ro','markersize',mSize);
plot(Mod_alphas(ex_set(2),mod_dt_ind),SU_ball_alphas(ex_set(2),direct_dt_ind,ball_ind),'ko','markersize',mSize);
line([0 1],[0 1]);
xlabel('Model-predicted alpha');
ylabel('Direct alpha');

%plot model fit for example neuron 1
npix = all_cell_data(ex_unit_ids(1),1).modFitParams.use_nPix_us;
flen = all_cell_data(ex_unit_ids(1),1).modFitParams.flen;
filt_mean = all_cell_data(ex_unit_ids(1),1).tune_props.avgRF_mean;
pix_dx = all_cell_data(ex_unit_ids(1),1).modFitParams.sp_dx;
pix_ax = (1:npix)*pix_dx; pix_ax = pix_ax - mean(pix_ax) - filt_mean;
dt = all_cell_data(ex_unit_ids(1),1).modFitParams.dt;
tax = (1:flen)*dt - dt/2;
figure(f2);clf
f2_props = plot_NMM_filters_1d(all_cell_data(ex_unit_ids(1),1).bestGQM,pix_ax,tax,[],f2,xr1,tr);

%plot model fit for example neuron 2
npix = all_cell_data(ex_unit_ids(2),1).modFitParams.use_nPix_us;
filt_mean = all_cell_data(ex_unit_ids(2),1).tune_props.avgRF_mean;
pix_dx = all_cell_data(ex_unit_ids(2),1).modFitParams.sp_dx;
pix_ax = (1:npix)*pix_dx; pix_ax = pix_ax - mean(pix_ax) - filt_mean;
flen = all_cell_data(ex_unit_ids(2),1).modFitParams.flen;
dt = all_cell_data(ex_unit_ids(2),1).modFitParams.dt;
tax = (1:flen)*dt - dt/2;

figure(f3);clf
f3_props = plot_NMM_filters_1d(all_cell_data(ex_unit_ids(2),1).bestGQM,pix_ax,tax,[],f3,xr2,tr);

%print out some properties of the stim tuning for these neurons
[RF_avg_width(ex_set(1)) RF_PRI(ex_set(1)) RF_ecc_avg(ex_set(1))]
[RF_avg_width(ex_set(2)) RF_PRI(ex_set(2)) RF_ecc_avg(ex_set(2))]


% fig_width = 4; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'mod_direct_alpha_withexamples.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

% fig_width = 2*f2_props.dims(2); rel_height = f2_props.dims(1)/f2_props.dims(2);
% figufy(f2);
% fname = [fig_dir sprintf('example_model_%d.pdf',SU_uset(ex_set(1)))];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% fig_width = 2*f3_props.dims(2); rel_height = f3_props.dims(1)/f3_props.dims(2);
% figufy(f3);
% fname = [fig_dir sprintf('example_model_%d.pdf',SU_uset(ex_set(2)))];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);




%% analyze validation based on simulated spiking
% close all
% min_mod_var = 1e-4; %minimum model-predicted total rate variance (otherwise, relative error blows up)
% modSim_uset = find(Mod_tot_vars(:,mod_dt_ind) > min_mod_var);
% fprintf('Found %d/%d with sufficient rate var\n',length(modSim_uset),length(SU_uset));
% 
% Mod_sim_alphas = nan(length(SU_uset),EP_params.sim_n_rpts);
% Mod_sim_ballvars = nan(length(SU_uset),EP_params.sim_n_rpts);
% for ii = 1:length(SU_uset)
%     Mod_sim_alphas(ii,:) = 1-arrayfun(@(x) x.eps_alphas(ball_ind),all_cell_data(SU_uset(ii),mod_dt_ind).mod_sim_stats);
%     Mod_sim_ballvars(ii,:) = arrayfun(@(x) x.eps_vars(ball_ind),all_cell_data(SU_uset(ii),mod_dt_ind).mod_sim_stats);
% end
% 
% %mean and SD of alpha and totvar estimates
% avg_sim_alphas = mean(Mod_sim_alphas,2);
% std_sim_alphas = std(Mod_sim_alphas,[],2);
% avg_sim_ballvars = mean(Mod_sim_ballvars,2);
% std_sim_ballvars = std(Mod_sim_ballvars,[],2);
% 
% %relative bias in alpha estimates 
% alpha_bias = (avg_sim_alphas - Mod_alphas(:,mod_dt_ind))./Mod_alphas(:,mod_dt_ind);
% 
% %relative uncertainty total variance estimates 
% bvar_rel_uncertainty = bsxfun(@rdivide,std_sim_ballvars,Mod_tot_vars(:,mod_dt_ind))*100; %relative uncertainty (%)
% alpha_rel_uncertainty = bsxfun(@rdivide,std_sim_alphas,Mod_alphas(:,mod_dt_ind))*100;
% 
% mod_dt = mod_bin_dts(mod_dt_ind);
% 
% f1 = figure();
% plot(Mod_tot_vars(modSim_uset,mod_dt_ind)/mod_dt^2,bvar_rel_uncertainty(modSim_uset),'.','markersize',mSize);
% suff_trials = modSim_uset(SU_nTrials(modSim_uset) > 50); %set of units with at least 50 trials
% hold on
% plot(Mod_tot_vars(suff_trials,mod_dt_ind)/mod_dt^2,bvar_rel_uncertainty(suff_trials),'r.','markersize',mSize);
% suff_trials = modSim_uset(SU_nTrials(modSim_uset) > 100); %set of units with at least 100 trials
% plot(Mod_tot_vars(suff_trials,mod_dt_ind)/mod_dt^2,bvar_rel_uncertainty(suff_trials),'k.','markersize',mSize);
% set(gca,'xscale','log');
% xlabel('Rate variance (Hz^2)');
% ylabel('Relative uncertainty (%)');
% legend('All units','At least 50 trials','At least 100 trials');
% 
% % fig_width = 4; rel_height = 0.8;
% % figufy(f1);
% % fname = [fig_dir 'Mod_rvar_unc.pdf'];
% % exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % % close(f1);
% % 
