clear all
close all

Expt_name = 'G087';
use_expt_num = 10;

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
% dir_prefix = '/media/NTlab_data1';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
cd(data_dir);


fname = sprintf('Expt%dClusterTimes.mat',use_expt_num);
load(fname);
for cc = 1:96
   if isfield(Clusters{cc},'marked')
       clust_mark(cc) = Clusters{cc}.marked;
   else
       clust_mark(cc) = 0;
   end
end


for pp = 1:96
    fprintf('Loading Expt%d P%d Spikes\n',use_expt_num,pp);
    fname = sprintf('Spikes/nby%s.p%dt%d.mat',Expt_name,pp,use_expt_num);
    load(fname);
    
    spk_values = double(Spikes.values);
    avg_wvfrm(:,pp) = mean(spk_values)*Spikes.maxv/Spikes.maxint;

end

%%
for pp = 1:96
    subplot(10,10,pp)
    if clust_mark(pp) == 2
    plot(avg_wvfrm(:,pp),'r')
    elseif clust_mark(pp) == 3
        plot(avg_wvfrm(:,pp),'g');
    elseif clust_mark(pp) == 4
        plot(avg_wvfrm(:,pp),'k')
    else
        plot(avg_wvfrm(:,pp));
    end
    ylim([-3 1]*1e-5);
    set(gca,'xticklabel',[],'yticklabel',[]);
    title(sprintf('Probe %d',pp));
end