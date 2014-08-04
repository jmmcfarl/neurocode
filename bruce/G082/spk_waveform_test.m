clear all
close all
cd /home/james/Data/bruce/G082

%%
% expt_num = 10;
% load(sprintf('Expt%dClusterTimesDetails.mat',expt_num));
% 
% %%
% % probe_num = 13;
% for probe_num = 1:96
%         
%     
%     %%
%     figure(1);
%     cur_inds = ClusterDetails{probe_num}.clst==2;
%     plot(ClusterDetails{probe_num}.xy(:,1),ClusterDetails{probe_num}.xy(:,2),'.','markersize',2)
%     hold on
%     plot(ClusterDetails{probe_num}.xy(cur_inds,1),ClusterDetails{probe_num}.xy(cur_inds,2),'r.','markersize',2)
%     yl = ylim();
%     line(ClusterDetails{probe_num}.crit([1 1]),yl,'color','k')
%     %%
%     load(sprintf('Expt%d.p%dFullV.mat',expt_num,probe_num));
%     dt =  FullV.samper;
%     Fs = 1/dt;
%     % Fs = 3e4;
%     V = double(FullV.V);
%     V = V*FullV.intscale(1)/FullV.intscale(2);
%     nparts = length(FullV.blklen);
%     dV = [];
%     %splice together multiple blocks
%     cur_pt = 1;
%     for pp = 1:nparts
%         cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
%         cur_range(cur_range > length(V)) = [];
%         curV = V(cur_range);
%         dV = [dV curV];
%         cur_pt = cur_pt + FullV.blklen(pp);
%     end
%     dV = dV(:);
%     
%     t_ax = [];
%     for pp = 1:nparts
%         cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
%         t_ax = [t_ax cur_t_ax];
%     end
%     t_ax(length(dV)+1:end) = [];
%     
% %     [filt_b,filt_a] = butter(2,100/(Fs/2),'high');
%     [filt_b,filt_a] = butter(2,[100 1e4]/(Fs/2));
%     V = filtfilt(filt_b,filt_a,dV);
%     %%
%     backlag = round(0.001*Fs);
%     forwardlag = round(0.002*Fs);
%     spkinds = ceil(interp1(t_ax,1:length(t_ax),ClusterDetails{probe_num}.t));
%     [avg_wv,lags] = get_event_trig_avg(V(:),spkinds,backlag,forwardlag);
%     [avg_wv_clust,lags] = get_event_trig_avg(V(:),spkinds(cur_inds),backlag,forwardlag);
%     figure(2);
%     plot(lags/Fs,avg_wv)
%     hold on
%     plot(lags/Fs,avg_wv_clust,'r')
%     
%     %%
%     ov_n_spks = length(spkinds);
%     clust_n_spks = sum(cur_inds);
%     fprintf('Rate > %.3f  Clustrate > %.3f\n',ov_n_spks/range(t_ax),clust_n_spks/range(t_ax));
%     
%     %%
%     pause
%     figure(1);clf;
%     figure(2);clf
% end
% 
% %%
% expt_num = 10;
% load(sprintf('Expt%dClusterTimesDetails.mat',expt_num));

%%
probe_num = 22;
for expt_num = 1:20
        
    load(sprintf('Expt%dClusterTimesDetails.mat',expt_num));
    load(sprintf('Expt%dClusterTimes.mat',expt_num));

    if expt_num==1
        base_ang = Clusters{probe_num}.angle;
        base_thresh = ClusterDetails{probe_num}.crit;
    end
    cur_rot = Clusters{probe_num}.angle - base_ang;
    rot_mat = [cos(cur_rot) -sin(cur_rot); sin(cur_rot) cos(cur_rot)];
    xy_data = ClusterDetails{probe_num}.xy*rot_mat;
%     cur_thresh = [ClusterDetails{probe_num}.crit 0]*rot_mat;
    %%
%     figure(1);
%     cur_inds = ClusterDetails{probe_num}.clst==2;
%     plot(ClusterDetails{probe_num}.xy(:,1),ClusterDetails{probe_num}.xy(:,2),'.','markersize',2)
%     hold on
%     plot(ClusterDetails{probe_num}.xy(cur_inds,1),ClusterDetails{probe_num}.xy(cur_inds,2),'r.','markersize',2)
%     yl = ylim();
%     line(ClusterDetails{probe_num}.crit([1 1]),yl,'color','k')
%     fprintf('Crit: %.3f\n',ClusterDetails{probe_num}.crit);
    figure(1);
    cur_inds = ClusterDetails{probe_num}.clst==2;
    plot(xy_data(:,1),xy_data(:,2),'.','markersize',2)
    hold on
    plot(xy_data(cur_inds,1),xy_data(cur_inds,2),'r.','markersize',2)
    yl = ylim();
    line(base_thresh([1 1]),yl,'color','k')
    fprintf('Crit: %.3f\n',ClusterDetails{probe_num}.crit);
    %%
    load(sprintf('Expt%d.p%dFullV.mat',expt_num,probe_num));
    dt =  FullV.samper;
    Fs = 1/dt;
    % Fs = 3e4;
    V = double(FullV.V);
    V = V*FullV.intscale(1)/FullV.intscale(2);
    nparts = length(FullV.blklen);
    dV = [];
    %splice together multiple blocks
    cur_pt = 1;
    for pp = 1:nparts
        cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
        cur_range(cur_range > length(V)) = [];
        curV = V(cur_range);
        dV = [dV curV];
        cur_pt = cur_pt + FullV.blklen(pp);
    end
    dV = dV(:);
    
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax cur_t_ax];
    end
    t_ax(length(dV)+1:end) = [];
    
%     [filt_b,filt_a] = butter(2,100/(Fs/2),'high');
    [filt_b,filt_a] = butter(2,[100 1e4]/(Fs/2));
    V = filtfilt(filt_b,filt_a,dV);
    %%
    backlag = round(0.001*Fs);
    forwardlag = round(0.002*Fs);
    spkinds = ceil(interp1(t_ax,1:length(t_ax),ClusterDetails{probe_num}.t));
    [avg_wv,lags] = get_event_trig_avg(V(:),spkinds,backlag,forwardlag);
    [avg_wv_clust,lags] = get_event_trig_avg(V(:),spkinds(cur_inds),backlag,forwardlag);
    figure(2);
    plot(lags/Fs,avg_wv)
    hold on
    plot(lags/Fs,avg_wv_clust,'r')
    
    %%
    ov_n_spks = length(spkinds);
    clust_n_spks = sum(cur_inds);
    fprintf('Rate > %.3f  Clustrate > %.3f\n',ov_n_spks/range(t_ax),clust_n_spks/range(t_ax));
    %%
    pause
    figure(1);clf;
    figure(2);clf
end

%%
for expt_num = 1:28
    expt_num
    load(sprintf('Expt%dClusterTimesDetails.mat',expt_num));
    load(sprintf('Expt%dClusterTimes.mat',expt_num));
for probe_num = 1:96;
        
    thresh(probe_num,expt_num) = ClusterDetails{probe_num}.crit;
       if isfield(Clusters{probe_num},'marked')
       clust_mark(probe_num,expt_num) = Clusters{probe_num}.marked;
       end
   clust_ang(probe_num,expt_num) = Clusters{probe_num}.angle;

end
end