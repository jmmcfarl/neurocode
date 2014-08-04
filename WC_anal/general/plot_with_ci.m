function [] = plot_with_ci(xdata,ymatrix,color)

ymean = nanmean(ymatrix);
numDat = size(ymatrix,1);
yup = ymean+2*nanstd(ymatrix)/sqrt(numDat);
ydown = ymean-2*nanstd(ymatrix)/sqrt(numDat);

plot(xdata,ymean,'Color',color,'linewidth',2)
hold on
plot(xdata,yup,'--','Color',color)
plot(xdata,ydown,'--','Color',color)