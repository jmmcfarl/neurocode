clear all
% close all

Expt_name = 'G081';
% dir_prefix = '/Volumes/james';
% dir_prefix = '~';
dir_prefix = '/media/NTlab_data1';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
cd(data_dir);
if ~exist('Spikes','dir')
    system('mkdir Spikes');
end

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
n_expts = length(Expts);
bad_expt = [];

Fs = 1/3.333307000208032e-05;
[filt_b,filt_a] = butter(4,100/(Fs/2),'high');
%%
for ee = 1:n_expts
    fprintf('Loading Expt%dClusterTimesDetails\n',ee);
    fname = sprintf('Expt%dClusterTimesDetails.mat',ee);
    load(fname);
    
    filename = sprintf('Expt%dFullVmean.mat',ee);
    load(filename);
    
    ll = 1;
    fprintf('Electrode %d of %d\n',ll,96);
    filename = sprintf('Expt%d.p%dFullV.mat',ee,ll);
    load(filename);
    V = double(FullV.V) + FullV.sumscale*sumv;
    
    nparts = length(FullV.blklen);
    cur_V = [];
    %splice together multiple blocks
    cur_pt = 1;
    for pp = 1:nparts
        cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
        cur_range(cur_range > length(V)) = [];
            tempV = filtfilt(filt_b,filt_a,V(cur_range));
            cur_V = [cur_V; tempV'];
        cur_pt = cur_pt + FullV.blklen(pp);
    end
%     cur_V = int16(cur_V);
    
    V_t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        V_t_ax = [V_t_ax; cur_t_ax'];
    end
    if length(V_t_ax) > size(cur_V,1)
        fprintf('%d extra time points, removing\n',length(V_t_ax) - size(cur_V,1));
        V_t_ax(size(cur_V,1)+1:end) = [];
    end
    
    for ll = 1:96
        fprintf('Electrode %d of %d\n',ll,96);
        filename = sprintf('Expt%d.p%dFullV.mat',ee,ll);
        load(filename);
        V = double(FullV.V) + FullV.sumscale*sumv;
%         V = double(FullV.V);
        
        nparts = length(FullV.blklen);
        cur_V = [];
        %splice together multiple blocks
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            tempV = filtfilt(filt_b,filt_a,V(cur_range));
%             tempV = V(cur_range);
            cur_V = [cur_V; tempV'];
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        cur_V = int16(cur_V);
        
        cur_spk_times = ClusterDetails{ll}.t;
        cur_spk_inds = round(interp1(V_t_ax,1:length(V_t_ax),cur_spk_times));
        time_err = abs(cur_spk_times - V_t_ax(cur_spk_inds)');
        bad_pts = find(time_err > 1/Fs);
        if ~isempty(bad_pts)
            fprintf('%d bad pts\n',length(bad_pts));
        end
        ev_mat = get_event_trig_mat(cur_V,cur_spk_inds,13,26);
        
        Spikes.values = int16(ev_mat);
        Spikes.maxv = FullV.intscale(1);
        Spikes.maxint = FullV.intscale(2);
        Spikes.times = V_t_ax(cur_spk_inds);
        spk_names = sprintf('Spikes/nby%s.p%dt%d.mat',Expt_name,ll,ee);
        fprintf('Saving %s\n',spk_names);
        save(spk_names,'Spikes');
        
    end
    
end


%%

