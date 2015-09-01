function [] = ebar_mat(x,M,ew,c)

n = size(M,1);
h = errorbar(x,squeeze(nanmean(M)),squeeze(nanstd(M))/sqrt(n),'color',c);
errorbar_tick(h,ew,'units');