function [] = boot_error_plot(Xmat,xax,nboot,color)

boot_avg = bootstrp(nboot,@nanmean,Xmat);
shadedErrorBar(xax,mean(boot_avg),std(boot_avg),{'color',color});