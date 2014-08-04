clear all
close all

Fs = 3e4;
dsf = 30;Fsd = Fs/dsf;
use_lfps = [1:16:96];
niqf = Fsd/2;
[bf,af] = butter(2,[60 120]/niqf);
%%
cd ~/Data/bruce/7_15_12/G034
all_V = [];
for ll = 1:length(use_lfps)
    fprintf('Electrode %d\n',ll);
    filename = sprintf('Expt%d.p%dFullV.mat',1,use_lfps(ll));
    load(filename);
    V = double(FullV.V);
%     V = decimate(V,dsf);
%     V = V*FullV.intscale(1)/FullV.intscale(2);
%     V = V(:);
    
    dV = [];
    cur_pt = 1;
        nparts = length(FullV.blklen);
    for pp = 1:nparts
        cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
        cur_range(cur_range > length(V)) = [];
        dV = [dV decimate(V(cur_range),dsf)];
        cur_pt = cur_pt + FullV.blklen(pp);
    end
    all_V(:,ll) = dV;
    all_Vf = filtfilt(bf,af,all_V);
end
t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
t_ax = downsample(t_ax,dsf);

%%
% params.Fs = Fsd;
% params.tapers = [6 11];
% win = 10;
window = round(Fsd*10);
clear Cxy
% for i = 2:length(use_lfps)
%     [C,phi,S12,S1,S2,f]=coherencysegc(all_V(:,1),all_V(:,i),win,params);
%     C_avg(i,:) = mean(C,2);
% end
for i = 2:length(use_lfps)
    %     [C,phi,S12,S1,S2,f]=coherencysegc(all_V2(:,1),all_V2(:,i),win,params);
    %     C_avg2(i,:) = mean(C,2);
    [Cxy(i,:),f] = mscohere(all_V(:,1),all_V(:,i),window,[],[],Fsd);
end
%%
use_expt = 6;
cd ~/Data/bruce/G075/
filename = sprintf('Expt%dFullVmean.mat',use_expt);
load(filename);

all_V2 = [];
for ll = 1:length(use_lfps)
    fprintf('Electrode %d\n',ll);
    filename = sprintf('Expt%d.p%dFullV.mat',use_expt,use_lfps(ll));
    load(filename);
    V = double(FullV.V);
    V = V + FullV.sumscale*sumv;
    %     V = decimate(V,dsf);
        V = V*FullV.intscale(1)/FullV.intscale(2);
    %     V = V(:);
    
    dV = [];
    cur_pt = 1;
    nparts = length(FullV.blklen);
    for pp = 1:nparts
        cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
        cur_range(cur_range > length(V)) = [];
        dV = [dV decimate(V(cur_range),dsf)];
        cur_pt = cur_pt + FullV.blklen(pp);
    end
    all_V2(:,ll) = dV;
    all_Vf2 = filtfilt(bf,af,all_V2);
end
    t_ax2 = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax2 = [t_ax2 downsample(cur_t_ax,dsf)];
    end
    t_ax2(size(all_V2,1)+1:end) = [];

%%
% params.Fs = Fsd;
% params.tapers = [6 11];
% win = 10;
window = round(Fsd*10);
clear Cxy2
for i = 2:length(use_lfps)
    %     [C,phi,S12,S1,S2,f]=coherencysegc(all_V2(:,1),all_V2(:,i),win,params);
    %     C_avg2(i,:) = mean(C,2);
    [Cxy2(i,:),f] = mscohere(all_V2(:,1),all_V2(:,i),window,[],[],Fsd);
end
% cur_seg = find(t_ax2 > 185 & t_ax2 < 188);
%     [Cxytemp,ftemp] = mscohere(all_V2(cur_seg,1),all_V2(cur_seg,2),[],[],[],Fsd);
