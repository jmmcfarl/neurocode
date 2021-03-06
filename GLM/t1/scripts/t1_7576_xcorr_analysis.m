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

%%
npairs  = ncells*(ncells-1)/2;
timeres = 0.0001; %in sec
tclen   = 0.05;
ncbins = tclen/timeres;
tmax    = max(max(cellfun(@(x)max(x),allspks)));
bin_ax = 0:timeres:tmax;
Nbins = length(bin_ax);
temp_ncbins = ncbins + 1;
t_ax = (-temp_ncbins:temp_ncbins)*timeres;
allccs = cell(ncells,ncells,nconds);
pairids = zeros(npairs,3);
for icell1 = 1:(ncells-1);
    icell1
    for icell2 = (icell1+1):ncells;
        icell2
        for n = 1:nconds
            %             cur_xcov = zeros(1,length(t_ax));
            %             for i = 1:nspks(n,icell1)
            %                 cur_xcov = cur_xcov + hist(allspks{n,icell2}-allspks{n,icell1}(i),t_ax);
            %             end
            %             cur_xcov([1 end]) = [];
            %             cur_xcov = (cur_xcov - nspks(n,icell1)*nspks(n,icell2)/Nbins)/sqrt((nspks(n,icell1) ...
            %                 - nspks(n,icell1)^2/Nbins)*(nspks(n,icell2) - nspks(n,icell2)^2/Nbins));
            %             allccs{icell1,icell2,n} = cur_xcov;
            
            strain1 = histc(allspks{n,icell1},bin_ax)';
            strain2 = histc(allspks{n,icell2},bin_ax);
            %             allccs{icell1,icell2,n} = xcov(strain1,strain2,ncbins,'biased'); %plot(allccs{ipair}); drawnow
            allccs{icell1,icell2,n} = xcov(strain1,strain2,ncbins,'coeff'); %plot(allccs{ipair}); drawnow
            
        end
    end
end;

% t_ax([1 end]) = [];

%%
temp_ncbins = ncbins + 1;
Nbins = length(bin_ax);
t_ax = (-temp_ncbins:temp_ncbins)*timeres;
cur_xcov = zeros(1,length(t_ax));
for i = 1:length(allspks{n,icell1})
    cur_xcov = cur_xcov + hist(allspks{n,icell2}-allspks{n,icell1}(i),t_ax);
end
cur_xcov([1 end]) = []; t_ax([1 end]) = [];
corrected_xcov = (cur_xcov - nspks(n,icell1)*nspks(n,icell2)/Nbins)/...
    sqrt((nspks(n,icell1) - nspks(n,icell1)^2/Nbins)*(nspks(n,icell2) - nspks(n,icell2)^2/Nbins));

%%
cd ~/Data/blanche/
save xcorrs7576_p1.mat allccs timeres tclen

%%
% taf        = rids{1}; %cellfun(@disp,cnames(tpn))
% twn        = rids{2};
% tpn        = rids{3}; %cellfun(@disp,cnames(tpn))
% tpt        = rids{4};
% tps        = rids{5}; %cellfun(@disp,cnames(tps))
% tns        = rids{6}; %cellfun(@disp,cnames(tns))
% 
% is_nat_contrast = ~cellfun(@(x) (x(end)=='5'),cnames);
% is_nat_mean = ~cellfun(@(x) (x(18)=='1'),cnames);
% stim_contrast = nan(size(is_nat_contrast));
% stim_contrast(~is_nat_contrast) = cellfun(@(x) str2num(x(end-2:end)),cnames(~is_nat_contrast));
% stim_mean = nan(size(is_nat_mean));
% stim_mean(~is_nat_mean) = cellfun(@(x) str2num(x(18:20)),cnames(~is_nat_mean));
% 
% clear used_conds
% 
% %DIFFERENT STIM GROUPS (POOLED ACROSS NS, PN, AND PS)
% % used_conds(1,:) = ismember(1:nconds,taf);
% % used_conds(2,:) = ismember(1:nconds,twn);
% % used_conds(3,:) = ismember(1:nconds,tpn);
% % used_conds(4,:) = ismember(1:nconds,tpt);
% % used_conds(5,:) = ismember(1:nconds,tps);
% % used_conds(6,:) = ~isnan(stim_contrast) & ~isnan(stim_mean) & ismember(1:nconds,tns);
% % used_conds(7,:) = ismember(1:nconds,[rids{3}(3:4) rids{4}(3:4) rids{5}(3:4) rids{6}(5:end)]);
% used_conds(1,:) = ismember(1:nconds,[rids{3}(3:4) rids{4}(3:4) rids{5}(3:4) rids{6}(5:end)]);
% 
% % %DIFFERENT CONTRAST GROUPS (POOLED ACROSS NS, PN, AND PS)
% % used_conds(1,:) = logical(stim_contrast == 15 & ismember(1:nconds,[tns tpn tps]));
% % used_conds(2,:) = logical(stim_contrast == 25 & ismember(1:nconds,[tns tpn tps]));
% % used_conds(3,:) = logical(stim_contrast == 35 & ismember(1:nconds,[tns tpn tps]));
% % used_conds(4,:) = logical(stim_contrast == 45 & ismember(1:nconds,[tns tpn tps]));
% % % used_conds(5,:) = logical(isnan(stim_contrast) & ismember(1:nconds,[tns tpn tps]));
% % used_conds(5,:) = logical(isnan(stim_contrast) & ~isnan(stim_mean) & ismember(1:nconds,[tns tpn tps]));
% 
% 
% % % % DIFFERENT MEAN GROUPS (POOLED ACROSS NS, PN, AND PS)
% % used_conds(1,:) = logical(stim_mean == 125 & ismember(1:nconds,[tns tpn tps]));
% % used_conds(2,:) = logical(isnan(stim_mean) & ismember(1:nconds,[tns tpn tps]));
% % % used_conds(2,:) = logical(isnan(stim_mean) & ~isnan(stim_contrast) & ismember(1:nconds,[tns tpn tps]));
% 
% n_cond_types = size(used_conds,1);
% cmap = colormap(jet (n_cond_types));
% 
% npairs  = ncells*(ncells-1)/2;
% ncbins = tclen/timeres;
% tmax    = max(max(cellfun(@(x)max(x),allspks)));
% bin_ax = 0:timeres:tmax;
% Nbins = length(bin_ax);
% t_ax = (-ncbins:ncbins)*timeres;
% pairids = zeros(npairs,3);
% 
% figure
% % cur_cell = 1;
% clear cur_xcov cur_rate
% for cur_cell = 1:26
%     cur_cell
%     cnt = 1;
%     for j = 1:ncells
%         if j ~= cur_cell
%             cur_xcov = zeros(size(used_conds,1),2*ncbins+1);
%             for m = 1:n_cond_types
%                 cur_cond_set = find(used_conds(m,:));
%                 for n = 1:length(cur_cond_set)
%                     if j > cur_cell
%                         cur_xcov(m,:) = cur_xcov(m,:) + allccs{cur_cell,j,cur_cond_set(n)};
%                     elseif j < cur_cell
%                         cur_xcov(m,:) = cur_xcov(m,:) + allccs{j,cur_cell,cur_cond_set(n)};
%                     end
%                 end
%                 cur_xcov(m,:) = cur_xcov(m,:)/length(cur_cond_set);
%                 cur_rate(m,cur_cell) = mean(nspks(cur_cond_set,cur_cell))/240;
%                 cur_rate(m,j) = mean(nspks(cur_cond_set,j))/240;
%                 subplot(5,5,cnt)
%                 %                                 title(sprintf('C1:%d  C2:%d',round(mean(nspks(cur_cond_set,cur_cell))),round(mean(nspks(cur_cond_set,j)))));
%                 if j > cur_cell
%                     used_xcov = cur_xcov(m,:);
%                 else
%                     used_xcov = fliplr(cur_xcov(m,:));
%                 end
%                 %                 plot(t_ax,used_xcov,'color',cmap(m,:),'linewidth',2); hold on
%                 plot(t_ax,used_xcov,'b'); hold on
%                 hold on
%                 smooth_cov = smooth(used_xcov,5,'lowess');
%                 plot(t_ax,smooth_cov,'r','linewidth',2)
%                 %                 plot(t_ax,
%                 axis tight; xlim([-0.005 0.005]);
%                 yl = ylim(); xl = xlim();
%                 line([0 0],yl,'color','k')
%                 line(xl,[0 0],'color','k')
%                 title(sprintf('C%d',j));
%             end
%             cnt = cnt + 1;
%         end
%     end
%     pause
%     clf
% end

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
% used_conds(1,:) = logical(stim_mean == 125 & ismember(1:nconds,[tns tpn tps]));
% used_conds(2,:) = logical(isnan(stim_mean) & ismember(1:nconds,[tns tpn tps]));
used_conds(1,:) = ismember(1:nconds,[rids{3}(3:4) rids{4}(3:4) rids{5}(3:4) rids{6}(5:end)]);

npairs  = ncells*(ncells-1)/2;
timeres = 0.0005; %in sec
tclen   = 0.015;
ncbins = tclen/timeres;
tmax    = max(max(cellfun(@(x)max(x),allspks)));
bin_ax = 0:timeres:tmax;
Nbins = length(bin_ax);
t_ax = (-ncbins:ncbins)*timeres;
cond_avg_ccs = cell(ncells,ncells);

%now compute precise xcorrs
htp_timeres = 0.0001; %in sec
htp_tclen   = 0.005;
htp_ncbins = tclen/htp_timeres;
htp_bin_ax = 0:htp_timeres:tmax;
htp_Nbins = length(htp_bin_ax);
htp_t_ax = (-htp_ncbins:htp_ncbins)*htp_timeres;
temp_t_ax = ((-htp_ncbins-1):(htp_ncbins+1))*htp_timeres;
htp_cond_avg_ccs = cell(ncells,ncells);

n_reps = 100;
jitter_window = 0.01;
cd ~/James_scripts/GLM/t1/rec7576_xcorrs/
for icell1 = 1:ncells-1
    for icell2 = (icell1+1):ncells
        fprintf('Pair %d-%d\n',icell1,icell2);
        for cond_t = 1:size(used_conds,1)
            cur_conds = find(used_conds(cond_t,:));
            cond_avg_ccs{icell1,icell2,cond_t} = zeros(1,length(t_ax));
            rand_ccs = zeros(length(cur_conds),n_reps,length(t_ax));
            for n = 1:length(cur_conds)
                fprintf('Cond %d of %d\n',n,length(cur_conds));
                cur_spks1 = allspks{cur_conds(n),icell1};
                cur_spks2 = allspks{cur_conds(n),icell2};
                
                %FASTER FOR BINNING > ~0.5ms
                strain1 = histc(cur_spks1,bin_ax);
                strain2 = histc(cur_spks2,bin_ax);
                cur_xcov = xcov(strain1,strain2,ncbins,'coeff'); %plot(allccs{ipair}); drawnow
                
                cond_avg_ccs{icell1,icell2,cond_t} = cond_avg_ccs{icell1,icell2,cond_t} + cur_xcov';
                
                for r = 1:n_reps
                    %                     rshift = ceil((rand-0.5)*jitter_window/timeres);
                    %                     strain2_shift = circshift(strain2,rshift);
                    
                    spk_shift = (rand(nspks(cur_conds(n),icell2),1)-0.5)*jitter_window;
                    strain2_shift = histc(cur_spks2+spk_shift,bin_ax);
                    rand_ccs(n,r,:) = xcov(strain1,strain2_shift,ncbins,'coeff');
                end
            end
            cond_avg_rand_ccs = squeeze(mean(rand_ccs));
            cond_avg_ccs{icell1,icell2,cond_t} = cond_avg_ccs{icell1,icell2,cond_t}/length(cur_conds);
            cond_rand_uci_ccs{icell1,icell2,cond_t} = prctile(cond_avg_rand_ccs,95);
            cond_rand_lci_ccs{icell1,icell2,cond_t} = prctile(cond_avg_rand_ccs,5);
            cond_rand_med_ccs{icell1,icell2,cond_t} = prctile(cond_avg_rand_ccs,50);
            cond_rand_avg_ccs{icell1,icell2,cond_t} = mean(cond_avg_rand_ccs);
            
            htp_cond_avg_ccs{icell1,icell2,cond_t} = zeros(1,length(htp_t_ax));
            for n = 1:length(cur_conds)
                cur_spks1 = allspks{cur_conds(n),icell1};
                cur_spks2 = allspks{cur_conds(n),icell2};
                
                %FASTER FOR BINNING < ~0.5ms
                cur_xcov = zeros(1,length(temp_t_ax));
                for i = 1:nspks(cur_conds(n),icell1)
                    cur_xcov = cur_xcov + hist(cur_spks2-cur_spks1(i),temp_t_ax);
                end
                cur_xcov([1 end]) = [];
                cur_xcov = (cur_xcov - nspks(cur_conds(n),icell1)*nspks(cur_conds(n),icell2)/...
                    htp_Nbins)/sqrt((nspks(cur_conds(n),icell1) - nspks(cur_conds(n),icell1)^2/htp_Nbins)...
                    *(nspks(cur_conds(n),icell2) - nspks(cur_conds(n),icell2)^2/htp_Nbins));
                
                htp_cond_avg_ccs{icell1,icell2,cond_t} = htp_cond_avg_ccs{icell1,icell2,cond_t} + cur_xcov;
            end
            htp_cond_avg_ccs{icell1,icell2,cond_t} = htp_cond_avg_ccs{icell1,icell2,cond_t}/length(cur_conds);
        end
        
            f2 = figure;
        set(f2,'PaperUnits','centimeters');
        set(f2, 'PaperSize', [30 55]);
        set(f2,'PaperPosition',[0,0,(get(f2,'PaperSize'))])
    subplot(2,1,1)
        set(gca,'fontname','arial','fontsize',14)
        plot(t_ax*1000,cond_avg_ccs{icell1,icell2,1},'b'); hold on
        hold on
        smooth_cov = smooth(cond_avg_ccs{icell1,icell2,1},7,'lowess');
        plot(t_ax*1000,smooth_cov,'r','linewidth',2)
        plot(t_ax*1000,cond_rand_med_ccs{icell1,icell2,1},'k')
        legend('Empirical','Smoothed','Jittered')
        plot(t_ax*1000,cond_rand_lci_ccs{icell1,icell2,1},'k--')
        plot(t_ax*1000,cond_rand_uci_ccs{icell1,icell2,1},'k--')
        axis tight; xlim([-0.01 0.01]*1000);
        yl = ylim(); xl = xlim();
        line([0 0],yl,'color','k')
        line(xl,[0 0],'color','k')
        xlabel('Time lag (ms)','fontsize',16)
        ylabel('Correlation coefficient','fontsize',16)
        title(sprintf('C%d x C%d',icell1,icell2),'fontsize',16);
        subplot(2,1,2)
        set(gca,'fontname','arial','fontsize',14)
        plot(htp_t_ax*1000,fliplr(htp_cond_avg_ccs{icell1,icell2,1}),'k','linewidth',2)
        xlim([-0.002 0.002]*1000)
        yl = ylim(); xl = xlim();
        line([0 0],yl,'color','k')
        line(xl,[0 0],'color','k')
        xlabel('Time lag (ms)','fontsize',16)
        ylabel('Correlation coefficient','fontsize',16)
         f_name = sprintf('C%d_C%d',icell1,icell2);
        print(f_name,'-dpng'); close all
       
    end
end

cd ~/Data/blanche/
save xcorrs7576_multires.mat allccs htp* timeres ncbins t_ax cond_*

%% COMPUTE AVG XCORRS ONLY FOR SPECIFIED CONDITIONS
cd ~/Data/blanche/
load xcorrs7576_multires

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
% used_conds(1,:) = logical(stim_mean == 125 & ismember(1:nconds,[tns tpn tps]));
% used_conds(2,:) = logical(isnan(stim_mean) & ismember(1:nconds,[tns tpn tps]));
used_conds(1,:) = ismember(1:nconds,[rids{3}(3:4) rids{4}(3:4) rids{5}(3:4) rids{6}(5:end)]);

tmax    = max(max(cellfun(@(x)max(x),allspks)));
ltp_timeres = 0.0025; %in sec
ltp_tclen   = 0.125;
ltp_ncbins = ltp_tclen/ltp_timeres;
ltp_bin_ax = 0:ltp_timeres:tmax;
ltp_Nbins = length(ltp_bin_ax);
ltp_t_ax = (-ltp_ncbins:ltp_ncbins)*ltp_timeres;
ltp_cond_avg_ccs = cell(ncells,ncells);
npairs  = ncells*(ncells-1)/2;

n_reps = 100;
jitter_window = 0.2;
cd ~/James_scripts/GLM/t1/rec7576_xcorrs/
for icell1 = 1:ncells-1
    for icell2 = (icell1+1):ncells
        fprintf('Pair %d-%d\n',icell1,icell2);
        for cond_t = 1:size(used_conds,1)
            cur_conds = find(used_conds(cond_t,:));
            ltp_cond_avg_ccs{icell1,icell2,cond_t} = zeros(1,length(ltp_t_ax));
            rand_ccs = zeros(length(cur_conds),n_reps,length(ltp_t_ax));
            for n = 1:length(cur_conds)
                fprintf('Cond %d of %d\n',n,length(cur_conds));
                cur_spks1 = allspks{cur_conds(n),icell1};
                cur_spks2 = allspks{cur_conds(n),icell2};
                
                %FASTER FOR BINNING > ~0.5ms
                strain1 = histc(cur_spks1,ltp_bin_ax);
                strain2 = histc(cur_spks2,ltp_bin_ax);
                cur_xcov = xcov(strain1,strain2,ltp_ncbins,'coeff'); %plot(allccs{ipair}); drawnow
                
                ltp_cond_avg_ccs{icell1,icell2,cond_t} = ltp_cond_avg_ccs{icell1,icell2,cond_t} + cur_xcov';
                
                for r = 1:n_reps
                    %                     rshift = ceil((rand-0.5)*jitter_window/timeres);
                    %                     strain2_shift = circshift(strain2,rshift);
                    
                    spk_shift = (rand(nspks(cur_conds(n),icell2),1)-0.5)*jitter_window;
                    strain2_shift = histc(cur_spks2+spk_shift,ltp_bin_ax);
                    rand_ccs(n,r,:) = xcov(strain1,strain2_shift,ltp_ncbins,'coeff');
                end
            end
            cond_avg_rand_ccs = squeeze(mean(rand_ccs));
            ltp_cond_avg_ccs{icell1,icell2,cond_t} = ltp_cond_avg_ccs{icell1,icell2,cond_t}/length(cur_conds);
            ltp_cond_rand_uci_ccs{icell1,icell2,cond_t} = prctile(cond_avg_rand_ccs,95);
            ltp_cond_rand_lci_ccs{icell1,icell2,cond_t} = prctile(cond_avg_rand_ccs,5);
            ltp_cond_rand_med_ccs{icell1,icell2,cond_t} = prctile(cond_avg_rand_ccs,50);
            ltp_cond_rand_avg_ccs{icell1,icell2,cond_t} = mean(cond_avg_rand_ccs);
            
        end
        
            f2 = figure('visible','off');
        set(f2,'PaperUnits','centimeters');
        set(f2, 'PaperSize', [30 55]);
        set(f2,'PaperPosition',[0,0,(get(f2,'PaperSize'))])
    subplot(3,1,1)
        set(gca,'fontname','arial','fontsize',14)
        plot(ltp_t_ax*1000,ltp_cond_avg_ccs{icell1,icell2,1},'b'); hold on
        hold on
        smooth_cov = smooth(ltp_cond_avg_ccs{icell1,icell2,1},7,'lowess');
        plot(ltp_t_ax*1000,smooth_cov,'r','linewidth',2)
        plot(ltp_t_ax*1000,ltp_cond_rand_avg_ccs{icell1,icell2,1},'k')
        legend('Empirical','Smoothed','Jittered')
        plot(ltp_t_ax*1000,ltp_cond_rand_lci_ccs{icell1,icell2,1},'k--')
        plot(ltp_t_ax*1000,ltp_cond_rand_uci_ccs{icell1,icell2,1},'k--')
        axis tight; xlim([-0.125 0.125]*1000);
        yl = ylim(); xl = xlim();
        line([0 0],yl,'color','k')
        line(xl,[0 0],'color','k')
        xlabel('Time lag (ms)','fontsize',16)
        ylabel('Correlation coefficient','fontsize',16)
        title(sprintf('C%d x C%d',icell1,icell2),'fontsize',16);
    subplot(3,1,2)
        set(gca,'fontname','arial','fontsize',14)
        plot(t_ax*1000,cond_avg_ccs{icell1,icell2,1},'b'); hold on
        hold on
        smooth_cov = smooth(cond_avg_ccs{icell1,icell2,1},7,'lowess');
        plot(t_ax*1000,smooth_cov,'r','linewidth',2)
        plot(t_ax*1000,cond_rand_avg_ccs{icell1,icell2,1},'k')
        legend('Empirical','Smoothed','Jittered')
        plot(t_ax*1000,cond_rand_lci_ccs{icell1,icell2,1},'k--')
        plot(t_ax*1000,cond_rand_uci_ccs{icell1,icell2,1},'k--')
        axis tight; xlim([-0.01 0.01]*1000);
        yl = ylim(); xl = xlim();
        line([0 0],yl,'color','k')
        line(xl,[0 0],'color','k')
        xlabel('Time lag (ms)','fontsize',16)
        ylabel('Correlation coefficient','fontsize',16)
        title(sprintf('C%d x C%d',icell1,icell2),'fontsize',16);
        subplot(3,1,3)
        set(gca,'fontname','arial','fontsize',14)
        plot(htp_t_ax*1000,fliplr(htp_cond_avg_ccs{icell1,icell2,1}),'k','linewidth',2)
        xlim([-0.002 0.002]*1000)
        yl = ylim(); xl = xlim();
        line([0 0],yl,'color','k')
        line(xl,[0 0],'color','k')
        xlabel('Time lag (ms)','fontsize',16)
        ylabel('Correlation coefficient','fontsize',16)
         f_name = sprintf('f2_C%d_C%d',icell1,icell2);
        print(f_name,'-dpng'); close all
       
    end
end

cd ~/Data/blanche/
save xcorrs7576_multires.mat allccs htp* timeres ncbins t_ax cond_*

%% MULTIRESOLUTION PLOTS
cd ~/Data/blanche/
load xcorrs7576_multires_meandep
ncells = 27;
figure;
for cur_cell = 1
    cur_cell
    cnt = 1;
    for j = cur_cell+1:ncells
        if j ~= cur_cell
            j
            subplot(2,1,1)
            set(gca,'fontname','arial','fontsize',14)
            plot(t_ax*1000,cond_avg_ccs{cur_cell,j,1},'b'); hold on
            hold on
            smooth_cov = smooth(cond_avg_ccs{cur_cell,j,1},7,'lowess');
            plot(t_ax*1000,smooth_cov,'r','linewidth',2)
            plot(t_ax*1000,cond_rand_med_ccs{cur_cell,j,1},'k')
            legend('Empirical','Smoothed','Jittered')
            plot(t_ax*1000,cond_rand_lci_ccs{cur_cell,j,1},'k--')
            plot(t_ax*1000,cond_rand_uci_ccs{cur_cell,j,1},'k--')
            axis tight; xlim([-0.015 0.015]*1000);
            yl = ylim(); xl = xlim();
            line([0 0],yl,'color','k')
            line(xl,[0 0],'color','k')
            xlabel('Time lag (ms)','fontsize',16)
            ylabel('Correlation coefficient','fontsize',16)
            title(sprintf('C%d x C%d',cur_cell,j),'fontsize',16);
            subplot(2,1,2)
            set(gca,'fontname','arial','fontsize',14)
            set(gca,'fontname','arial','fontsize',14)
            plot(t_ax*1000,cond_avg_ccs{cur_cell,j,2},'b'); hold on
            hold on
            smooth_cov = smooth(cond_avg_ccs{cur_cell,j,2},7,'lowess');
            plot(t_ax*1000,smooth_cov,'r','linewidth',2)
            plot(t_ax*1000,cond_rand_med_ccs{cur_cell,j,2},'k')
            legend('Empirical','Smoothed','Jittered')
            plot(t_ax*1000,cond_rand_lci_ccs{cur_cell,j,2},'k--')
            plot(t_ax*1000,cond_rand_uci_ccs{cur_cell,j,2},'k--')
              axis tight; xlim([-0.015 0.015]*1000);
            yl = ylim(); xl = xlim();
            line([0 0],yl,'color','k')
            line(xl,[0 0],'color','k')
            xlabel('Time lag (ms)','fontsize',16)
            ylabel('Correlation coefficient','fontsize',16)
            pause
            clf
        end
        cnt = cnt + 1;
    end
    %     pause
    %     clf
end

%% MULTIRESOLUTION PLOTS
cd ~/Data/blanche/
load xcorrs7576_multires
ncells = 27;
figure;
% cur_cell = 10;
% j = 14;
for cur_cell = 5
    cur_cell
    cnt = 1;
    for j = cur_cell+1:ncells
        if j ~= cur_cell
            j
            subplot(2,1,1)
            set(gca,'fontname','arial','fontsize',14)
            plot(t_ax*1000,cond_avg_ccs{cur_cell,j},'b'); hold on
            hold on
            smooth_cov = smooth(cond_avg_ccs{cur_cell,j},7,'lowess');
            plot(t_ax*1000,smooth_cov,'r','linewidth',2)
            plot(t_ax*1000,cond_rand_med_ccs{cur_cell,j},'k')
            legend('Empirical','Smoothed','Jittered')
            plot(t_ax*1000,cond_rand_lci_ccs{cur_cell,j},'k--')
            plot(t_ax*1000,cond_rand_uci_ccs{cur_cell,j},'k--')
            axis tight; xlim([-0.015 0.015]*1000);
            yl = ylim(); xl = xlim();
            line([0 0],yl,'color','k')
            line(xl,[0 0],'color','k')
            xlabel('Time lag (ms)','fontsize',16)
            ylabel('Correlation coefficient','fontsize',16)
            title(sprintf('C%d x C%d',cur_cell,j),'fontsize',16);
            subplot(2,1,2)
            set(gca,'fontname','arial','fontsize',14)
            plot(htp_t_ax*1000,fliplr(htp_cond_avg_ccs{cur_cell,j}),'k','linewidth',2)
            xlim([-0.002 0.002]*1000)
            yl = ylim(); xl = xlim();
            line([0 0],yl,'color','k')
            line(xl,[0 0],'color','k')
            xlabel('Time lag (ms)','fontsize',16)
            ylabel('Correlation coefficient','fontsize',16)
            pause
            clf
        end
        cnt = cnt + 1;
    end
    %     pause
    %     clf
end




%% MULTIRESOLUTION MUTLICELL PLOTS
cd ~/Data/blanche/
load xcorrs7576_multires
c1 = 1;
c2 = [2];

figure;
for j = 1:length(c2)
    subplot(length(c2),2,2*(j-1)+1)
    set(gca,'fontname','arial','fontsize',14)
    if c2(j) > c1
        plot(t_ax*1000,cond_avg_ccs{c1,c2(j)},'b'); hold on
        hold on
        smooth_cov = smooth(cond_avg_ccs{c1,c2(j)},7,'lowess');
        plot(t_ax*1000,smooth_cov,'r','linewidth',2)
        plot(t_ax*1000,cond_rand_med_ccs{c1,c2(j)},'k')
%         legend('Empirical','Smoothed','Jittered')
        plot(t_ax*1000,cond_rand_lci_ccs{c1,c2(j)},'k--')
        plot(t_ax*1000,cond_rand_uci_ccs{c1,c2(j)},'k--')
    else
        plot(t_ax*1000,fliplr(cond_avg_ccs{c2(j),c1}),'b'); hold on
        hold on
        smooth_cov = fliplr(smooth(cond_avg_ccs{c2(j),c1},7,'lowess'));
        plot(t_ax*1000,smooth_cov,'r','linewidth',2)
        plot(t_ax*1000,fliplr(cond_rand_med_ccs{c2(j),c1}),'k')
%         legend('Empirical','Smoothed','Jittered')
        plot(t_ax*1000,fliplr(cond_rand_lci_ccs{c2(j),c1}),'k--')
        plot(t_ax*1000,fliplr(cond_rand_uci_ccs{c2(j),c1}),'k--')        
    end
    axis tight; xlim([-0.010 0.010]*1000);
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time lag (ms)','fontsize',16)
    ylabel('Correlation coefficient','fontsize',16)
    title(sprintf('C%d x C%d',c1,c2(j)),'fontsize',16);
    
    subplot(length(c2),2,2*(j-1)+2)
    set(gca,'fontname','arial','fontsize',14)
    if c2(j) > c1
    plot(htp_t_ax*1000,fliplr(htp_cond_avg_ccs{c1,c2(j)}),'k','linewidth',2)
    else
     plot(htp_t_ax*1000,htp_cond_avg_ccs{c2(j),c1},'k','linewidth',2)       
    end
    axis tight
    xlim([-0.002 0.002]*1000)
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time lag (ms)','fontsize',16)
    ylabel('Correlation coefficient','fontsize',16)
end
cd ~/James_scripts/GLM/t1/
%% Get xcorr mat for just the data used for model fitting
xcorr_mat = cell(26,26);
cur_cond_set = find(used_conds(7,:));

t_ax = (-ncbins:ncbins)*timeres;

for i = 1:25
    for j = i+1:26
        cur_xcov = zeros(1,2*ncbins+1);
        for n = 1:length(cur_cond_set)
            cur_xcov = cur_xcov + allccs{i,j,cur_cond_set(n)};
        end
        xcorr_mat{i,j} = cur_xcov/length(cur_cond_set);
    end
end
cd '/Users/James/James_scripts/GLM/t1/'
save xcorr_mat_7576_usedconds used_conds xcorr_mat t_ax timeres ncbins

%%
figure
c1 = 1;
c2 = 25;
% c1 = 1;
% c2 = 25;
% for c2 = 2:26;
cur_xcov = zeros(size(used_conds,1),2*ncbins+1);
for m = 1:n_cond_types
    cur_cond_set = find(used_conds(m,:));
    for n = 1:length(cur_cond_set)
        cur_xcov(m,:) = cur_xcov(m,:) + allccs{c1,c2,cur_cond_set(n)};
    end
    cur_xcov(m,:) = cur_xcov(m,:)/length(cur_cond_set);
    subplot(2,1,1)
    plot(t_ax,cur_xcov(m,:),'color',cmap(m,:),'linewidth',2); hold on
    subplot(2,1,2)
    plot(1,cur_rate(m,c1),'o','color',cmap(m,:),'linewidth',2); hold on
    plot(2,cur_rate(m,c2),'o','color',cmap(m,:),'linewidth',2); hold on
end
subplot(2,1,1)
xlim([-0.005 0.005])
xlabel('Time lag (s)','fontsize',14)
ylabel('Correlation','fontsize',14)
legend('Mn-125','Mn-NAT')
subplot(2,1,2)
xlim([0.5 2.5])
xlabel('Cell number','fontsize',14)
ylabel('Firing rate (Hz)','fontsize',14)
yl = ylim;
yl(1) = 0;
ylim(yl)

%     pause
%     clf
% end
%%
n_cond_types = size(used_conds,1);
figure
cmap = colormap(jet(n_cond_types));
% c1 = 12;
% c2 = [10 14 19];
% c1 = 3;
% c2 = [1 4 7];
% c1 = 2;
% c2 = [16 17 18 20];
c1 = 1;
c2 = [3 4 25];
nc2 = length(c2);
for cc = 1:nc2
    cur_xcov = zeros(size(used_conds,1),2*ncbins+1);
    for m = 1:n_cond_types
        cur_cond_set = find(used_conds(m,:));
        for n = 1:length(cur_cond_set)
            if c2(cc) > c1
                cur_xcov(m,:) = cur_xcov(m,:) + allccs{c1,c2(cc),cur_cond_set(n)};
            else
                cur_xcov(m,:) = cur_xcov(m,:) + fliplr(allccs{c2(cc),c1,cur_cond_set(n)});
            end
        end
        cur_xcov(m,:) = cur_xcov(m,:)/length(cur_cond_set);
        subplot(nc2,1,cc)
        plot(t_ax,cur_xcov(m,:),'color',cmap(m,:),'linewidth',2); hold on
        xlim([-0.005 0.005])
    end
end
subplot(nc2,1,1)
legend('C-15','C-25','C-35','C-45','C-NAT')
% legend('M-125','M-NAT')

%%
% c1 = 1;
% c2 = 25;
% bin_ax = 0:timeres:tmax;
% cc_confs = nan(nconds,2,2*ncbins+1);
% nreps = 500;
% for n = 1:nconds
%     strain1 = histc(allspks{n,c1},bin_ax)';
%     temp_conf = nan(nreps,2*ncbins+1);
%     for r = 1:nreps
%         rshift = rand*120;
%         strain2 = histc(mod(allspks{n,c2}+rshift,120),bin_ax);
%         temp_conf(r,:) = xcov(strain1,strain2,ncbins,'coeff');
%     end
%     cc_confs(n,1,:) = prctile(temp_conf,5);
%     cc_confs(n,2,:) = prctile(temp_conf,95);
% end
%

