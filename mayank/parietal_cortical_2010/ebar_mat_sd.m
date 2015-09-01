function [] = ebar_mat_sd(x,M,ew,c)

n = size(M,1);
h = errorbar(x,squeeze(nanmean(M)),squeeze(nanstd(M)),'color',c);
errorbar_tick(h,ew,'units');