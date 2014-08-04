function new_autoclust = try_mult_gmms(probe_num,stored_spkxy,stored_times,autoclust,...
    poss_n_comps,n_reps,Expt_name,target_probes,target_expts)

n_expts = size(autoclust,1);
new_autoclust = autoclust;

for ee = target_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    fname = sprintf('Spikes/nby%s.p%dt%d.mat',Expt_name,target_probes(probe_num),ee);
    load(fname);
    
    n_clst_spks = length(stored_times{ee,probe_num});
    if n_clst_spks ~= length(Spikes.times)
        use_set = find(ismember(Spikes.times,stored_times{ee,probe_num}));
        use_set2 = find(ismember(stored_times{ee,probe_num},Spikes.times));
        fprintf('Warning, no match for %d of %d spikes\n',length(setdiff(1:n_clst_spks,use_set)),n_clst_spks);
    else
        use_set = 1:n_clst_spks;
        use_set2 = 1:n_clst_spks;
    end
    spk_values = double(Spikes.values(use_set,:));
    
    for nc = 1:length(poss_n_comps)
        n_comps = poss_n_comps(nc);
        clust_gmm = gmdistribution.fit(stored_spkxy{ee,probe_num},n_comps,'Replicates',n_reps,'Regularize',1e-12,'Options',statset('MaxIter',1000));
        [clust_idx,nlogl] = cluster(clust_gmm,stored_spkxy{ee,probe_num});
        
        cur_avg_wvfrm = zeros(40,n_comps);
        for nn = 1:n_comps
            cur_avg_wvfrm(:,nn) = mean(spk_values(clust_idx(use_set2)==nn,:));
        end
        
        smallest_vals = min(cur_avg_wvfrm);
        [~,clust_ord] = sort(smallest_vals,'ascend');
        temp_idx = clust_idx;
        temp_idx(clust_idx == clust_ord(1)) = 1;
        temp_idx(clust_idx ~= clust_ord(1)) = 2;
        clust_idx = temp_idx;
        
        used_clust_mean = clust_gmm.mu(clust_ord(1),:);
        used_clust_cov = squeeze(clust_gmm.Sigma(:,:,clust_ord(1)));
        other_clust_mean = mean(clust_gmm.mu(clust_ord(2:end),:),1);
        other_clust_cov = mean(clust_gmm.Sigma(:,:,clust_ord(2:end)),3);
        
        mean_diff = used_clust_mean - other_clust_mean;
        D1 = mean_diff*inv(other_clust_cov)*mean_diff';
        D2 = mean_diff*inv(used_clust_cov)*mean_diff';
        mahal_d = sqrt(2./((1./D1)+(1./D2)));
        
        %         D1 = mahal(stored_spkxy{ee,probe_num}(clust_idx == 1,:),stored_spkxy{ee,probe_num}(clust_idx == 2,:));
        %         D2 = mahal(stored_spkxy{ee,probe_num}(clust_idx == 2,:),stored_spkxy{ee,probe_num}(clust_idx == 1,:));
        %         mahal_d = sqrt(2./((1./mean(D1))+(1./mean(D2))));
        
        mv1 = mean(stored_spkxy{ee,probe_num}(clust_idx == 1,:));
        mv2 = mean(stored_spkxy{ee,probe_num}(clust_idx == 2,:));
        cov1 = cov(stored_spkxy{ee,probe_num}(clust_idx == 1,:));
        cov2 = cov(stored_spkxy{ee,probe_num}(clust_idx == 2,:));
        D1 = (mv1-mv2)*inv(cov2)*(mv1-mv2)';
        D2 = (mv1-mv2)*inv(cov1)*(mv1-mv2)';
        tmahal_d = sqrt(2./((1./D1)+(1./D2)));
        
        avg_wvfrm(:,1) = mean(spk_values(clust_idx(use_set2) == 1,:));
        avg_wvfrm(:,2) = mean(spk_values(clust_idx(use_set2) == 2,:));
        std_wvfrm(:,1) = std(spk_values(clust_idx(use_set2) == 1,:));
        std_wvfrm(:,2) = std(spk_values(clust_idx(use_set2) == 2,:));
        avg_wvfrm = avg_wvfrm*Spikes.maxv/Spikes.maxint;
        std_wvfrm = std_wvfrm*Spikes.maxv/Spikes.maxint;
        
        tempclust(ee,nc).gmm = clust_gmm;
        tempclust(ee,nc).n_comps = n_comps;
        tempclust(ee,nc).idx = clust_idx;
        tempclust(ee,nc).avg_wvfrm = avg_wvfrm;
        tempclust(ee,nc).std_wvfrm = std_wvfrm;
        tempclust(ee,nc).nlogl = nlogl;
        tempclust(ee,nc).mahal_d = mahal_d;
        tempclust(ee,nc).tmahal_d = tmahal_d;
    end
end
%% NOW EXAMINE RESULTS AND PICK BEST FITS IN EACH CASE
n_cols = ceil(sqrt(length(poss_n_comps)));
n_rows = ceil(length(poss_n_comps)/n_cols);
wrat = n_cols/n_rows;

f1 = figure;
f1_pos = get(f1,'Position');
f1_pos(3) = f1_pos(4)*1.25*wrat;
set(f1,'Position',f1_pos);
f2 = figure;
f2_pos = get(f2,'Position');
f2_pos(3) = f2_pos(4)*1.25*wrat;
f2_pos(2) = f2_pos(2) - 500;
set(f2,'Position',f2_pos);
for ee = target_expts
    figure(f1);clf
    figure(f2);clf
    for nc = 1:length(poss_n_comps)
        figure(f1);
        subplot(n_rows,n_cols,nc)
        plot(stored_spkxy{ee,probe_num}(tempclust(ee,nc).idx == 1,1),stored_spkxy{ee,probe_num}(tempclust(ee,nc).idx == 1,2),'r.','markersize',0.5);
        hold on
        plot(stored_spkxy{ee,probe_num}(tempclust(ee,nc).idx == 2,1),stored_spkxy{ee,probe_num}(tempclust(ee,nc).idx == 2,2),'k.','markersize',0.5);
        set(gca,'xticklabel',[],'yticklabel',[]);
        title(sprintf('Expt %d, %d comps',ee,poss_n_comps(nc)));
        axis tight
        title(sprintf('%d comp Mah = %.3f',poss_n_comps(nc),tempclust(ee,nc).tmahal_d));
        
        figure(f2);
        subplot(n_rows,n_cols,nc)
        [handles, details] = DensityPlot_jmm(stored_spkxy{ee,probe_num}(:,1),stored_spkxy{ee,probe_num}(:,2),'sqrtsc');
        hold on
        z = pdf(tempclust(ee,nc).gmm,[details.x(:) details.y(:)]);
        z = sqrt(z);
        contour(details.x, details.y, reshape(z,size(details.x)),10,'w','linewidth',0.5);
        set(gca,'xticklabel',[],'yticklabel',[],'ydir','normal');
        title(sprintf('Expt %d, %d comps',ee,poss_n_comps(nc)));
        
    end
    
    use_linesep = 0; use_robust = 0; refit = 0;
    have_char = 0;
    while have_char == 0
        ui = input(sprintf('Expt %d use which model? (q to exit, 0 use old, l# line-sep, r# repeat robust)\n',ee),'s');
        if ~isempty(ui) 
            have_char = 1;
        end
        if isempty(ui) && length(poss_n_comps) == 1
           have_char = 1;
           ui = num2str(poss_n_comps); %assume the only option was entered
        end
        if have_char == 0
            fprintf('Invalid input\n');
        end
    end
    if strcmpi(ui,'q')
        break
    end
    if strcmpi(ui(1),'l')
        use_linesep = 1;
        refit = 1;
        if length(poss_n_comps) > 1
            ui = ui(2);
        else
            ui = num2str(poss_n_comps);
        end
    elseif strcmpi(ui(1),'r')
        refit = 1;
        use_robust = 1;
        if length(poss_n_comps) > 1
            ui = ui(2);
        else
            ui = num2str(poss_n_comps);
        end
    end
    ui = str2num(ui);
    if ui == 0 %use old clust
        new_autoclust(ee,probe_num).man_code = 0;
    elseif refit == 0
        while ~ismember(ui,poss_n_comps)
            fprintf('%d not a possible choice, re-enter\n',ui);
            ui = input(sprintf('Expt %d use which model? (#only)\n',ee),'s');
            ui = str2num(ui);
        end
        fprintf('Setting to %d comp-GMM\n',ui);
        cur_clust = tempclust(ee,poss_n_comps == ui);
        new_autoclust(ee,probe_num).man_code = 1;
        new_autoclust(ee,probe_num).gmm = cur_clust.gmm;
        new_autoclust(ee,probe_num).idx = cur_clust.idx;
        new_autoclust(ee,probe_num).avg_wvfrm = cur_clust.avg_wvfrm;
        new_autoclust(ee,probe_num).std_wvfrm = cur_clust.std_wvfrm;
        new_autoclust(ee,probe_num).nlogl = cur_clust.nlogl;
        new_autoclust(ee,probe_num).mahal_d = cur_clust.mahal_d;
        new_autoclust(ee,probe_num).tmahal_d = cur_clust.tmahal_d;
        new_autoclust(ee,probe_num).dprime = nan;
        new_autoclust(ee,probe_num).n_comps = cur_clust.n_comps;
    else
            while refit == 1
                cur_n_comps = ui;
                if use_robust == 1
                    newtempclust = robust_clust(Expt_name,target_probes(probe_num),ee,stored_times{ee,probe_num},stored_spkxy{ee,probe_num},cur_n_comps);
                elseif use_linesep == 1
                    newtempclust = linesep_clust(Expt_name,target_probes(probe_num),ee,stored_times{ee,probe_num},stored_spkxy{ee,probe_num},cur_n_comps);
                end
                figure(f1);
                subplot(n_rows,n_cols,find(poss_n_comps == cur_n_comps))
                plot(stored_spkxy{ee,probe_num}(newtempclust.idx == 1,1),stored_spkxy{ee,probe_num}(newtempclust.idx == 1,2),'r.','markersize',0.5);
                hold on
                plot(stored_spkxy{ee,probe_num}(newtempclust.idx == 2,1),stored_spkxy{ee,probe_num}(newtempclust.idx == 2,2),'k.','markersize',0.5);
                set(gca,'xticklabel',[],'yticklabel',[]);
                title(sprintf('Expt %d, %d comps',ee,poss_n_comps(nc)));
                axis tight
                title(sprintf('%d comp Mah = %.3f',poss_n_comps(nc),newtempclust.tmahal_d));
                figure(f2);
                subplot(n_rows,n_cols,find(poss_n_comps == cur_n_comps))
                [handles, details] = DensityPlot_jmm(stored_spkxy{ee,probe_num}(:,1),stored_spkxy{ee,probe_num}(:,2),'sqrtsc');
                hold on
                z = pdf(newtempclust.gmm,[details.x(:) details.y(:)]);
                z = sqrt(z);
                contour(details.x, details.y, reshape(z,size(details.x)),10,'w','linewidth',0.5);
                set(gca,'xticklabel',[],'yticklabel',[],'ydir','normal');
                title(sprintf('Expt %d, %d comps',ee,cur_n_comps));
                
                ni = input('Accept? (y/# for yes, r to retry robust, l to retry linesep, a to abort)\n','s');
                if strcmpi(ni,'a')
                    refit = 0;
                    new_autoclust(ee,probe_num).man_code = 0;
                elseif strcmpi(ni,'y') || any(strcmp({'1' '2' '3' '4' '5' '6' '7' '8' '9'},ni)) || isempty(ni)
                    fprintf('Using this %d-comp model\n',cur_n_comps);
                    refit = 0;
                    cur_clust = newtempclust;
                    new_autoclust(ee,probe_num).man_code = 1;
                    new_autoclust(ee,probe_num).gmm = cur_clust.gmm;
                    new_autoclust(ee,probe_num).idx = cur_clust.idx;
                    new_autoclust(ee,probe_num).avg_wvfrm = cur_clust.avg_wvfrm;
                    new_autoclust(ee,probe_num).std_wvfrm = cur_clust.std_wvfrm;
                    new_autoclust(ee,probe_num).nlogl = cur_clust.nlogl;
                    new_autoclust(ee,probe_num).mahal_d = cur_clust.mahal_d;
                    new_autoclust(ee,probe_num).tmahal_d = cur_clust.tmahal_d;
                    new_autoclust(ee,probe_num).dprime = nan;
                    new_autoclust(ee,probe_num).n_comps = cur_clust.n_comps;
                elseif strcmpi(ni,'r')
                    use_robust = 1;
                    use_linesep = 0;
                elseif strcmpi(ni,'l')
                    use_robust = 0;
                    use_linesep = 1;
                else
                    fprintf('Invalid input\n');
                    ti = input('Abort or retry (a/r)?\n','s');
                    if strcmpi(ti,'a')
                        refit = 0;
                        new_autoclust(ee,probe_num).man_code = 0;
                    end
                end
        end
    end
end
end

%%
function tempclust = robust_clust(ename,probe,expt,stored_times,stored_xy,n_comps)
n_reps = 10;
mah_thresh = 7;

fname = sprintf('Spikes/nby%s.p%dt%d.mat',ename,probe,expt);
load(fname);

n_clst_spks = length(stored_times);
if n_clst_spks ~= length(Spikes.times)
    use_set = find(ismember(Spikes.times,stored_times));
    use_set2 = find(ismember(stored_times,Spikes.times));
    fprintf('Warning, no match for %d of %d spikes\n',length(setdiff(1:n_clst_spks,use_set)),n_clst_spks);
else
    use_set = 1:n_clst_spks;
    use_set2 = 1:n_clst_spks;
end
spk_values = double(Spikes.values(use_set,:));


clust_gmm = gmdistribution.fit(stored_xy,n_comps,'Replicates',n_reps,'Regularize',1e-12,'Options',statset('MaxIter',1000));
[clust_idx,nlogl] = cluster(clust_gmm,stored_xy);
cur_mah = mahal(clust_gmm,stored_xy);
near_mah = min(cur_mah,[],2);
not_outliers = find(near_mah < mah_thresh);
% km_idx = kmeans(stored_xy(not_outliers),n_comps,'distance','cityblock','Replicates',100);
% clust_gmm = gmdistribution.fit(stored_xy(not_outliers,:),n_comps,'Start',km_idx,'Regularize',1e-12,'Options',statset('MaxIter',1000));
% [clust_idx,nlogl] = cluster(clust_gmm,stored_xy);
clust_gmm = gmdistribution.fit(stored_xy(not_outliers,:),n_comps,'Replicates',n_reps,'Regularize',1e-12,'Options',statset('MaxIter',1000));
[clust_idx,nlogl] = cluster(clust_gmm,stored_xy);

cur_avg_wvfrm = zeros(40,n_comps);
for nn = 1:n_comps
    cur_avg_wvfrm(:,nn) = mean(spk_values(clust_idx(use_set2)==nn,:));
end

smallest_vals = min(cur_avg_wvfrm);
[~,clust_ord] = sort(smallest_vals,'ascend');
temp_idx = clust_idx;
temp_idx(clust_idx == clust_ord(1)) = 1;
temp_idx(clust_idx ~= clust_ord(1)) = 2;
clust_idx = temp_idx;

used_clust_mean = clust_gmm.mu(clust_ord(1),:);
used_clust_cov = squeeze(clust_gmm.Sigma(:,:,clust_ord(1)));
other_clust_mean = mean(clust_gmm.mu(clust_ord(2:end),:),1);
other_clust_cov = mean(clust_gmm.Sigma(:,:,clust_ord(2:end)),3);

mean_diff = used_clust_mean - other_clust_mean;
D1 = mean_diff*inv(other_clust_cov)*mean_diff';
D2 = mean_diff*inv(used_clust_cov)*mean_diff';
mahal_d = sqrt(2./((1./D1)+(1./D2)));

mv1 = mean(stored_xy(clust_idx == 1,:));
mv2 = mean(stored_xy(clust_idx == 2,:));
cov1 = cov(stored_xy(clust_idx == 1,:));
cov2 = cov(stored_xy(clust_idx == 2,:));
D1 = (mv1-mv2)*inv(cov2)*(mv1-mv2)';
D2 = (mv1-mv2)*inv(cov1)*(mv1-mv2)';
tmahal_d = sqrt(2./((1./D1)+(1./D2)));

avg_wvfrm(:,1) = mean(spk_values(clust_idx(use_set2) == 1,:));
avg_wvfrm(:,2) = mean(spk_values(clust_idx(use_set2) == 2,:));
std_wvfrm(:,1) = std(spk_values(clust_idx(use_set2) == 1,:));
std_wvfrm(:,2) = std(spk_values(clust_idx(use_set2) == 1,:));
avg_wvfrm = avg_wvfrm*Spikes.maxv/Spikes.maxint;
std_wvfrm = std_wvfrm*Spikes.maxv/Spikes.maxint;

tempclust.gmm = clust_gmm;
tempclust.n_comps = n_comps;
tempclust.idx = clust_idx;
tempclust.avg_wvfrm = avg_wvfrm;
tempclust.std_wvfrm = std_wvfrm;
tempclust.nlogl = nlogl;
tempclust.mahal_d = mahal_d;
tempclust.tmahal_d = tmahal_d;
end

%%
function tempclust = linesep_clust(ename,probe,expt,stored_times,stored_xy,n_comps)

fname = sprintf('Spikes/nby%s.p%dt%d.mat',ename,probe,expt);
load(fname);

n_clst_spks = length(stored_times);
if n_clst_spks ~= length(Spikes.times)
    use_set = find(ismember(Spikes.times,stored_times));
    use_set2 = find(ismember(stored_times,Spikes.times));
    fprintf('Warning, no match for %d of %d spikes\n',length(setdiff(1:n_clst_spks,use_set)),n_clst_spks);
else
    use_set = 1:n_clst_spks;
    use_set2 = 1:n_clst_spks;
end
spk_values = double(Spikes.values(use_set,:));

acc = 0;

ftemp = figure;
while acc == 0
    plot(stored_xy(:,1),stored_xy(:,2),'k.','markersize',0.5);
    fprintf('Create line with two points\n');
    [x,y] = ginput(2);
    sep_vec = [x(2)-x(1); y(2)-y(1)];
    orth_vec = [0 1;-1 0]*sep_vec;
    proj_data = stored_xy*orth_vec;
    proj_thresh = [x(1) y(1)]*orth_vec;
    above_thresh = find(proj_data > proj_thresh);
    hold on
    plot(stored_xy(above_thresh,1),stored_xy(above_thresh,2),'r.','markersize',0.5);
    %     ui = input('Accept (y)\n','s');
    %     if strcmpi(ui,'y')
    %         acc = 1;
    %     end
    pause(0.25);
    acc = 1;
end
close(ftemp);
init_idx = ones(size(stored_xy,1),1);
init_idx(above_thresh) = 2;

clust_gmm = gmdistribution.fit(stored_xy,n_comps,'Start',init_idx,'Regularize',1e-12,'Options',statset('MaxIter',1000));
[clust_idx,nlogl] = cluster(clust_gmm,stored_xy);

cur_avg_wvfrm = zeros(40,n_comps);
for nn = 1:n_comps
    cur_avg_wvfrm(:,nn) = mean(spk_values(clust_idx(use_set2)==nn,:));
end

smallest_vals = min(cur_avg_wvfrm);
[~,clust_ord] = sort(smallest_vals,'ascend');
temp_idx = clust_idx;
temp_idx(clust_idx == clust_ord(1)) = 1;
temp_idx(clust_idx ~= clust_ord(1)) = 2;
clust_idx = temp_idx;

used_clust_mean = clust_gmm.mu(clust_ord(1),:);
used_clust_cov = squeeze(clust_gmm.Sigma(:,:,clust_ord(1)));
other_clust_mean = mean(clust_gmm.mu(clust_ord(2:end),:),1);
other_clust_cov = mean(clust_gmm.Sigma(:,:,clust_ord(2:end)),3);

mean_diff = used_clust_mean - other_clust_mean;
D1 = mean_diff*inv(other_clust_cov)*mean_diff';
D2 = mean_diff*inv(used_clust_cov)*mean_diff';
mahal_d = sqrt(2./((1./D1)+(1./D2)));

mv1 = mean(stored_xy(clust_idx == 1,:));
mv2 = mean(stored_xy(clust_idx == 2,:));
cov1 = cov(stored_xy(clust_idx == 1,:));
cov2 = cov(stored_xy(clust_idx == 2,:));
D1 = (mv1-mv2)*inv(cov2)*(mv1-mv2)';
D2 = (mv1-mv2)*inv(cov1)*(mv1-mv2)';
tmahal_d = sqrt(2./((1./D1)+(1./D2)));

avg_wvfrm(:,1) = mean(spk_values(clust_idx(use_set2) == 1,:));
avg_wvfrm(:,2) = mean(spk_values(clust_idx(use_set2) == 2,:));
std_wvfrm(:,1) = std(spk_values(clust_idx(use_set2) == 1,:));
std_wvfrm(:,2) = std(spk_values(clust_idx(use_set2) == 1,:));
avg_wvfrm = avg_wvfrm*Spikes.maxv/Spikes.maxint;
std_wvfrm = std_wvfrm*Spikes.maxv/Spikes.maxint;

tempclust.gmm = clust_gmm;
tempclust.n_comps = n_comps;
tempclust.idx = clust_idx;
tempclust.avg_wvfrm = avg_wvfrm;
tempclust.std_wvfrm = std_wvfrm;
tempclust.nlogl = nlogl;
tempclust.mahal_d = mahal_d;
tempclust.tmahal_d = tmahal_d;
end
