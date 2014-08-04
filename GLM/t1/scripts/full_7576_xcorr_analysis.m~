clear all;
cd ~/Data/blanche/rec_75/matlabdata/
stf75=load('stimfiles75.mat');
load spksegsRec75.mat;
load stdparsRec75;
cd ~/Data/blanche/rec_76/matlabdata/
stf76=load('stimfileNames76.mat');
load spksegsRec76.mat;

allspks  = [spksegs75,spksegs76]';

%merge cell 10 and 14 (really same cell) to create cell 27
cell10_spks = allspks(:,10);
cell14_spks = allspks(:,14);
temp_comb = cell(size(allspks,1),1);
for i = 1:size(allspks,1)
    temp_comb{i} = unique([cell10_spks{i};cell14_spks{i}]);
end
allspks = [allspks temp_comb];

nspks = cellfun(@length,allspks);

[nconds,ncells]   = size(nspks);

cnames   = [stf75.stimfiles,stf76.stimfiles];
cellfun(@disp,cnames)

keys = {'af','wn','pn','pt','ps','ns'}
rids = cellfun(@(tkey)...
    find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
    keys,'UniformOutput',0);

%% COMPUTE AVG XCORRS ONLY FOR SPECIFIED CONDITIONS
is_nat_contrast = ~cellfun(@(x) (x(end)=='5'),cnames);
is_nat_mean = ~cellfun(@(x) (x(18)=='1'),cnames);
stim_contrast = nan(size(is_nat_contrast));
stim_contrast(~is_nat_contrast) = cellfun(@(x) str2num(x(end-2:end)),cnames(~is_nat_contrast));
stim_mean = nan(size(is_nat_mean));
stim_mean(~is_nat_mean) = cellfun(@(x) str2num(x(18:20)),cnames(~is_nat_mean));


%DIFFERENT STIM GROUPS (POOLED ACROSS NS, PN, AND PS)
clear used_conds
% used_conds = find(ismember(1:nconds,[rids{3}(3:4) rids{4}(3:4) rids{5}(3:4) rids{6}(5:end)]));
taf        = rids{1}; %cellfun(@disp,cnames(tpn))
twn        = rids{2};
tpn        = rids{3}; %cellfun(@disp,cnames(tpn))
tpt        = rids{4};
tps        = rids{5}; %cellfun(@disp,cnames(tps))
tns        = rids{6}; %cellfun(@disp,cnames(tns))
all_conds = 1:nconds;
used_conds = find(ismember(all_conds,[rids{3}(3:4) rids{4}(3:4) rids{5}(3:4) rids{6}(3:end)]));
con_vals = [15 25 35 45 nan];

tmax    = max(max(cellfun(@(x)max(x),allspks)));
npairs  = ncells*(ncells-1)/2;

ltp_timeres = 0.0025; %in sec
ltp_tclen   = 0.100;
ltp_ncbins = ltp_tclen/ltp_timeres;
ltp_bin_ax = 0:ltp_timeres:tmax;
ltp_Nbins = length(ltp_bin_ax);
ltp_t_ax = (-ltp_ncbins:ltp_ncbins)*ltp_timeres;
ltp_cond_avg_ccs = cell(ncells,ncells);
ltp_jitter_window = 0.1;

mtp_timeres = 0.0005; %in sec
mtp_tclen   = 0.01;
mtp_ncbins = mtp_tclen/mtp_timeres;
mtp_bin_ax = 0:mtp_timeres:tmax;
mtp_Nbins = length(mtp_bin_ax);
mtp_t_ax = (-mtp_ncbins:mtp_ncbins)*mtp_timeres;
mtp_cond_avg_ccs = cell(ncells,ncells);
mtp_jitter_window = 0.01;

htp_timeres = 0.0001; %in sec
htp_tclen   = 0.005;
htp_ncbins = htp_tclen/htp_timeres;
htp_bin_ax = 0:htp_timeres:tmax;
htp_Nbins = length(htp_bin_ax);
htp_t_ax = (-htp_ncbins:htp_ncbins)*htp_timeres;
temp_t_ax = ((-htp_ncbins-1):(htp_ncbins+1))*htp_timeres;
htp_cond_avg_ccs = cell(ncells,ncells);

n_reps = 100;
for icell1 = 1:ncells-1
    for icell2 = (icell1+1):ncells
        fprintf('Pair %d-%d\n',icell1,icell2);
        
        
        ltp_cond_avg_ccs{icell1,icell2} = zeros(1,length(ltp_t_ax));
        ltp_stim_dep_ccs = zeros(length(all_conds),length(ltp_t_ax));
        rand_ccs = zeros(length(all_conds),n_reps,length(ltp_t_ax));
        for n = 1:length(all_conds)
            fprintf('Cond %d of %d\n',n,length(all_conds));
            cur_spks1 = allspks{all_conds(n),icell1};
            cur_spks2 = allspks{all_conds(n),icell2};
            
            strain1 = histc(cur_spks1,ltp_bin_ax);
            strain2 = histc(cur_spks2,ltp_bin_ax);
            cur_xcov = xcov(strain1,strain2,ltp_ncbins,'coeff'); %plot(allccs{ipair}); drawnow
            
            ltp_cond_avg_ccs{icell1,icell2} = ltp_cond_avg_ccs{icell1,icell2} + cur_xcov';
            ltp_stim_dep_ccs(n,:) = cur_xcov;
            for r = 1:n_reps
                spk_shift = (rand(nspks(all_conds(n),icell2),1)-0.5)*ltp_jitter_window;
                strain2_shift = histc(cur_spks2+spk_shift,ltp_bin_ax);
                rand_ccs(n,r,:) = xcov(strain1,strain2_shift,ltp_ncbins,'coeff');
            end
        end
        cond_avg_rand_ccs = squeeze(mean(rand_ccs));
        ltp_cond_avg_ccs{icell1,icell2} = ltp_cond_avg_ccs{icell1,icell2}/length(all_conds);
        ltp_cond_rand_uci_ccs{icell1,icell2} = prctile(cond_avg_rand_ccs,95);
        ltp_cond_rand_lci_ccs{icell1,icell2} = prctile(cond_avg_rand_ccs,5);
        ltp_cond_rand_med_ccs{icell1,icell2} = prctile(cond_avg_rand_ccs,50);
        ltp_cond_rand_avg_ccs{icell1,icell2} = mean(cond_avg_rand_ccs);
        ltp_cond_zavg_ccs{icell1,icell2} = (ltp_cond_avg_ccs{icell1,icell2} - mean(cond_avg_rand_ccs))./std(cond_avg_rand_ccs);
        
        mtp_cond_avg_ccs{icell1,icell2} = zeros(1,length(mtp_t_ax));
        mtp_stim_dep_ccs = zeros(length(all_conds),length(mtp_t_ax));
        rand_ccs = zeros(length(all_conds),n_reps,length(mtp_t_ax));
        for n = 1:length(all_conds)
            fprintf('Cond %d of %d\n',n,length(all_conds));
            cur_spks1 = allspks{all_conds(n),icell1};
            cur_spks2 = allspks{all_conds(n),icell2};
            
            strain1 = histc(cur_spks1,mtp_bin_ax);
            strain2 = histc(cur_spks2,mtp_bin_ax);
            cur_xcov = xcov(strain1,strain2,mtp_ncbins,'coeff'); %plot(allccs{ipair}); drawnow
            
            mtp_cond_avg_ccs{icell1,icell2} = mtp_cond_avg_ccs{icell1,icell2} + cur_xcov';
            mtp_stim_dep_ccs(n,:) = cur_xcov;
            for r = 1:n_reps
                spk_shift = (rand(nspks(all_conds(n),icell2),1)-0.5)*mtp_jitter_window;
                strain2_shift = histc(cur_spks2+spk_shift,mtp_bin_ax);
                rand_ccs(n,r,:) = xcov(strain1,strain2_shift,mtp_ncbins,'coeff');
            end
        end
        cond_avg_rand_ccs = squeeze(mean(rand_ccs));
        mtp_cond_avg_ccs{icell1,icell2} = mtp_cond_avg_ccs{icell1,icell2}/length(all_conds);
        mtp_cond_rand_uci_ccs{icell1,icell2} = prctile(cond_avg_rand_ccs,95);
        mtp_cond_rand_lci_ccs{icell1,icell2} = prctile(cond_avg_rand_ccs,5);
        mtp_cond_rand_med_ccs{icell1,icell2} = prctile(cond_avg_rand_ccs,50);
        mtp_cond_rand_avg_ccs{icell1,icell2} = mean(cond_avg_rand_ccs);
        mtp_cond_zavg_ccs{icell1,icell2} = (mtp_cond_avg_ccs{icell1,icell2} - mean(cond_avg_rand_ccs))./std(cond_avg_rand_ccs);
        
        
        htp_cond_avg_ccs{icell1,icell2} = zeros(1,length(htp_t_ax));
        htp_stim_dep_ccs = zeros(length(all_conds),length(htp_t_ax));
        for n = 1:length(all_conds)
            cur_spks1 = allspks{all_conds(n),icell1};
            cur_spks2 = allspks{all_conds(n),icell2};
            
            cur_xcov = zeros(1,length(temp_t_ax));
            for i = 1:nspks(all_conds(n),icell1)
                cur_xcov = cur_xcov + hist(cur_spks2-cur_spks1(i),temp_t_ax);
            end
            cur_xcov([1 end]) = [];
            cur_xcov = (cur_xcov - nspks(all_conds(n),icell1)*nspks(all_conds(n),icell2)/...
                htp_Nbins)/sqrt((nspks(all_conds(n),icell1) - nspks(all_conds(n),icell1)^2/htp_Nbins)...
                *(nspks(all_conds(n),icell2) - nspks(all_conds(n),icell2)^2/htp_Nbins));
            
            htp_cond_avg_ccs{icell1,icell2} = htp_cond_avg_ccs{icell1,icell2} + fliplr(cur_xcov);
            htp_stim_dep_ccs(n,:) = fliplr(cur_xcov);
        end
        htp_cond_avg_ccs{icell1,icell2} = htp_cond_avg_ccs{icell1,icell2}/length(all_conds);
        
        
        cd ~/James_scripts/GLM/t1/rec7576_xcorrs/
        f2 = figure('visible','off');
        set(f2,'PaperUnits','centimeters');
        set(f2, 'PaperSize', [55 45]);
        set(f2,'PaperPosition',[0,0,(get(f2,'PaperSize'))])
        subplot(3,3,1)
        set(gca,'fontname','arial','fontsize',14)
        plot(ltp_t_ax*1000,ltp_cond_avg_ccs{icell1,icell2},'b'); hold on
        hold on
        smooth_cov = smooth(ltp_cond_avg_ccs{icell1,icell2},7,'lowess');
        plot(ltp_t_ax*1000,smooth_cov,'r','linewidth',2)
        plot(ltp_t_ax*1000,ltp_cond_rand_avg_ccs{icell1,icell2},'k')
        legend('Empirical','Smoothed','Jittered')
        plot(ltp_t_ax*1000,ltp_cond_rand_lci_ccs{icell1,icell2},'k--')
        plot(ltp_t_ax*1000,ltp_cond_rand_uci_ccs{icell1,icell2},'k--')
        axis tight; xlim([-0.1 0.1]*1000);
        yl = ylim(); xl = xlim();
        line([0 0],yl,'color','k')
        line(xl,[0 0],'color','k')
        xlabel('Time lag (ms)','fontsize',16)
        ylabel('Correlation coefficient','fontsize',16)
        title(sprintf('C%d x C%d',icell1,icell2),'fontsize',16);
        subplot(3,3,4)
        set(gca,'fontname','arial','fontsize',14)
        plot(mtp_t_ax*1000,mtp_cond_avg_ccs{icell1,icell2},'b'); hold on
        hold on
        smooth_cov = smooth(mtp_cond_avg_ccs{icell1,icell2},7,'lowess');
        plot(mtp_t_ax*1000,smooth_cov,'r','linewidth',2)
        plot(mtp_t_ax*1000,mtp_cond_rand_avg_ccs{icell1,icell2},'k')
        legend('Empirical','Smoothed','Jittered')
        plot(mtp_t_ax*1000,mtp_cond_rand_lci_ccs{icell1,icell2},'k--')
        plot(mtp_t_ax*1000,mtp_cond_rand_uci_ccs{icell1,icell2},'k--')
        axis tight; xlim([-0.01 0.01]*1000);
        yl = ylim(); xl = xlim();
        line([0 0],yl,'color','k')
        line(xl,[0 0],'color','k')
        xlabel('Time lag (ms)','fontsize',16)
        ylabel('Correlation coefficient','fontsize',16)
        title(sprintf('C%d x C%d',icell1,icell2),'fontsize',16);
        subplot(3,3,7)
        set(gca,'fontname','arial','fontsize',14)
        plot(htp_t_ax*1000,htp_cond_avg_ccs{icell1,icell2},'k','linewidth',2)
        axis tight; xlim([-0.002 0.002]*1000)
        yl = ylim(); xl = xlim();
        line([0 0],yl,'color','k')
        line(xl,[0 0],'color','k')
        xlabel('Time lag (ms)','fontsize',16)
        ylabel('Correlation coefficient','fontsize',16)
        
        cmap = colormap(jet(length(con_vals)));
        for i = 1:length(con_vals)-1
            cur_conds = find(stim_contrast==con_vals(i));
            subplot(3,3,2)
            plot(ltp_t_ax*1000,mean(ltp_stim_dep_ccs(cur_conds,:)),'color',cmap(i,:),'linewidth',1)
            hold on
            subplot(3,3,5)
            plot(mtp_t_ax*1000,mean(mtp_stim_dep_ccs(cur_conds,:)),'color',cmap(i,:),'linewidth',1)
            hold on
            subplot(3,3,8)
            plot(htp_t_ax*1000,mean(htp_stim_dep_ccs(cur_conds,:)),'color',cmap(i,:),'linewidth',1);
            hold on
        end
        cur_conds = find(isnan(stim_contrast));
        subplot(3,3,2)
        plot(ltp_t_ax*1000,mean(ltp_stim_dep_ccs(cur_conds,:)),'color',cmap(end,:),'linewidth',1)
        hold on
        subplot(3,3,5)
        plot(mtp_t_ax*1000,mean(mtp_stim_dep_ccs(cur_conds,:)),'color',cmap(end,:),'linewidth',1)
        hold on
        subplot(3,3,8)
        plot(htp_t_ax*1000,mean(htp_stim_dep_ccs(cur_conds,:)),'color',cmap(end,:),'linewidth',1);
        hold on
        subplot(3,3,2)
        xlim([-0.1 0.1]*1000);
        plot(ltp_t_ax*1000,ltp_cond_avg_ccs{icell1,icell2},'color','k','linewidth',2)
        subplot(3,3,5)
        xlim([-0.01 0.01]*1000);
        plot(mtp_t_ax*1000,mtp_cond_avg_ccs{icell1,icell2},'color','k','linewidth',2)
        subplot(3,3,8)
        xlim([-0.002 0.002]*1000);
        plot(htp_t_ax*1000,htp_cond_avg_ccs{icell1,icell2},'color','k','linewidth',2)
        
        
        cmap = colormap(jet(length(rids)));
        for i = 1:length(rids)
            cur_conds = rids{i};
            subplot(3,3,3)
            plot(ltp_t_ax*1000,mean(ltp_stim_dep_ccs(cur_conds,:)),'color',cmap(i,:),'linewidth',1)
            hold on
            subplot(3,3,6)
            plot(mtp_t_ax*1000,mean(mtp_stim_dep_ccs(cur_conds,:)),'color',cmap(i,:),'linewidth',1)
            hold on
            subplot(3,3,9)
            plot(htp_t_ax*1000,mean(htp_stim_dep_ccs(cur_conds,:)),'color',cmap(i,:),'linewidth',1);
            hold on
        end
        subplot(3,3,3)
        xlim([-0.1 0.1]*1000);
        plot(ltp_t_ax*1000,ltp_cond_avg_ccs{icell1,icell2},'color','k','linewidth',2)
        subplot(3,3,6)
        xlim([-0.01 0.01]*1000);
        plot(mtp_t_ax*1000,mtp_cond_avg_ccs{icell1,icell2},'color','k','linewidth',2)
        subplot(3,3,9)
        xlim([-0.002 0.002]*1000);
        plot(htp_t_ax*1000,htp_cond_avg_ccs{icell1,icell2},'color','k','linewidth',2)
        
        
        f_name = sprintf('ff_C%d_C%d',icell1,icell2);
        print(f_name,'-dpng'); close all
        cd ~/Data/blanche/
        save xcorrs7576_fullmultires.mat htp* mtp* ltp*
        
    end
end

%%
cnt = 1;
ltp_zero_bin = find(ltp_t_ax==0);
no_zero_bin = 1:length(ltp_t_ax); no_zero_bin(ltp_zero_bin) = [];
for icell1 = 1:ncells-1
    for icell2 = (icell1+1):ncells
        pair_id(cnt,1) = icell1;
        pair_id(cnt,2) = icell2;
        
        cur_ltp_corr = ltp_cond_zavg_ccs{icell1,icell2};
        cur_ltp_corr = interp1(no_zero_bin,cur_ltp_corr(no_zero_bin),1:length(ltp_t_ax));
        [a,b] = max(abs(cur_ltp_corr));
        ltp_max_zcorr(cnt) = cur_ltp_corr(b);
        ltp_max_zlag(cnt) = ltp_t_ax(b)*1000;
        ltp_max_corr(cnt) = ltp_cond_avg_ccs{icell1,icell2}(b);
        
         [a,b] = max(abs(mtp_cond_zavg_ccs{icell1,icell2}));
        mtp_max_zcorr(cnt) = mtp_cond_zavg_ccs{icell1,icell2}(b);
        mtp_max_zlag(cnt) = mtp_t_ax(b)*1000;
        mtp_max_corr(cnt) = mtp_cond_avg_ccs{icell1,icell2}(b);
       
        cur_htp_corr = htp_cond_avg_ccs{icell1,icell2};
        [a,b] = max(abs(cur_htp_corr));
        htp_max_corr(cnt) = cur_htp_corr(b);
        htp_max_lag(cnt) = htp_t_ax(b)*1000;
        htp_max_zcorr(cnt) = (cur_htp_corr(b)-mean(cur_htp_corr))/std(cur_htp_corr);
        
        cnt = cnt + 1;
    end
end

sig_ltp_corrs = find(abs(ltp_max_zcorr) > 5);

%locate artifact peaks
thresh_corr = 6e-3;
min_width = 4;
pot_art = find(htp_max_corr > thresh_corr);
peak_width = nan(size(pot_art));
for i = 1:length(pot_art)
    cur_pairs = pair_id(pot_art(i),:);   
    cur_corr = htp_cond_avg_ccs{cur_pairs(1),cur_pairs(2)};
    [a,b] = max(cur_corr);
    r_w = find(cur_corr(b:end) < htp_max_corr(pot_art(i))*0.25,1,'first')-1;
    r_l = b-find(cur_corr(1:b) < htp_max_corr(pot_art(i))*0.25,1,'last');
    peak_width(i) = r_w + r_l;
end
artifacts = pot_art(peak_width < min_width);

n_pairs = (ncells-1)*ncells/2;
used_pairs = setdiff(1:n_pairs,artifacts);

%%
figure
plot(ltp_max_zcorr(used_pairs),mtp_max_zcorr(used_pairs),'.')

figure
plot(ltp_max_corr(used_pairs),mtp_max_corr(used_pairs),'.')

figure
plot(ltp_max_corr(used_pairs),htp_max_corr(used_pairs),'.')

figure
plot(ltp_max_zcorr(used_pairs),htp_max_zcorr(used_pairs),'.')

figure
plot(mtp_max_zcorr(used_pairs),htp_max_zcorr(used_pairs),'.')

%%
zero_lag = find(mtp_max_zlag(used_pairs) == 0);
nzero_lag = find(abs(mtp_max_zlag(used_pairs)) >= 0.5 & abs(mtp_max_zlag(used_pairs)) <= 4);
