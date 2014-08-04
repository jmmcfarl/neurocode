function [] = dataViewer_run_hist2(dataw,data8,spktimes,spkAmp)

dsf = 8;
Fsd = 2016/8;
segSize = 30;
t_axis = [1:length(data8)]/Fsd;

spktimes = round(spktimes/dsf);
spkAmp = zscore(spkAmp);
spkhist = hist(spktimes,[1:length(data8)]);
spkrate = jmm_smooth_1d(spkhist,20)/Fsd;
isis = [0;diff(spktimes)];
isis(isis>100) = 100;
subplot(3,1,1)
plot(t_axis,data8,'linewidth',2)
subplot(3,1,2)
plot(t_axis,dataw,'r','linewidth',2)
subplot(3,1,3)
plot(t_axis(spktimes),isis,'k.-')
% 
% plot(t_axis,int_thresh_8,'k','linewidth',2)
% plot(t_axis(up_times8),int_thresh_8(up_times8),'go')
% plot(t_axis(down_times8),int_thresh_8(down_times8),'ro')

% plot(t_axis(spktimes),ones(size(spktimes))*2,'k.-')
spkrate = zscore(spkrate);
% plot(t_axis,spkrate,'c')

numSegs = floor(max(t_axis)/segSize);

for i = 1:numSegs
    subplot(3,1,1)
    xlim([segSize*(i-1)+1 segSize*i])
    ylim([-3 3])
   subplot(3,1,2)
       xlim([segSize*(i-1)+1 segSize*i])
    ylim([-3 3])

   subplot(3,1,3)
    xlim([segSize*(i-1)+1 segSize*i])

   pause
    
    
end