function [] = dataViewer_run_hist(data,data8,int_thresh,int_thresh_8,up_times,up_times8,down_times,down_times8)

Fsd = 2016/8;
segSize = 50;
t_axis = [1:length(data)]/Fsd;

subplot(2,1,1)
plot(t_axis,data,'linewidth',2)
hold on
plot(t_axis,int_thresh,'k','linewidth',2)

plot(t_axis(up_times),int_thresh(up_times),'go')
plot(t_axis(down_times),int_thresh(down_times),'ro')

subplot(2,1,2)
plot(t_axis,data8,'linewidth',2)
hold on
plot(t_axis,int_thresh_8,'k','linewidth',2)
plot(t_axis(up_times8),int_thresh_8(up_times8),'go')
plot(t_axis(down_times8),int_thresh_8(down_times8),'ro')

numSegs = floor(max(t_axis)/segSize);

for i = 1:numSegs
    subplot(2,1,1)
    xlim([segSize*(i-1)+1 segSize*i])
    ylim([-3 3])
    subplot(2,1,2)
    xlim([segSize*(i-1)+1 segSize*i])
    ylim([-3 3])
    pause
    
end