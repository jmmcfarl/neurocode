function [] = matrix_errorbar_plot(Xmat,xax,color)

n_trials = size(Xmat,1);
avg = squeeze(nanmean(Xmat));
sem = squeeze(nanstd(Xmat))/sqrt(n_trials);

shadedErrorBar(xax,avg,sem,{'color',color});
