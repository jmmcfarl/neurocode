close all
clear all
fig_dir = '/home/james/Analysis/bruce/saccade_modulation/';

%% LOAD JBE HORIZONTAL
sname = 'sacStimProc_sta';
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
rmfield_list = {};

all_hor_data = [];
all_hor_tdata = [];
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(save_dir)
    load(sname)
    
    ucells = arrayfun(@(x) length(x.ModData),sacStimProc) > 0;
    cur_data = sacStimProc(ucells);
    cur_data = rmfields(cur_data,rmfield_list);
    [cur_data.expt_num] = deal(Expt_num);
    [cur_data.bar_ori] = deal(0);
    
    all_hor_data = cat(2,all_hor_data,cur_data);
    cur_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,cur_data);

end

%%
xx = linspace(-5,5,500);
pp = normpdf(xx,0,1);
pp = pp/sum(pp);

alpha = 1;
beta_range = linspace(0.01,1.5,500);
theta_range = linspace(-6,3,200);
% theta_range = 0;

info = nan(length(beta_range),length(theta_range));
for bb = 1:length(beta_range)
    for tt = 1:length(theta_range)
        rr = alpha*log(1+exp(beta_range(bb)*xx+theta_range(tt)));
        mr = sum(pp.*rr);
        info(bb,tt) = sum(pp.*rr.*log2(rr/mr))/mr;
        mrate(bb,tt) = mr;
    end
end

%%

betas = [all_hor_data(:).ov_beta];
thetas = [all_hor_data(:).ov_theta];

figure;
pcolor(beta_range,theta_range,sqrt(info)');shading flat
hold on
plot(betas,thetas,'wo','linewidth',2);
plot(mean(betas),mean(thetas),'g*','linewidth',2);
plot(1,mean(thetas),'r*','linewidth',2);
plot(0.1,mean(thetas),'m*','linewidth',2);

figure;
rr = alpha*log(1+exp(mean(betas)*xx+mean(thetas)));

plot(xx,rr);
yy = ylim();
hold on
plot(xx,pp/max(pp)*yy(2)*0.8,'k');

rr = alpha*log(1+exp(1*xx+mean(thetas)));
plot(xx,rr/max(rr)*yy(2),'r');
rr = alpha*log(1+exp(0.1*xx+mean(thetas)));
rr = rr-min(rr);
plot(xx,rr/max(rr)*yy(2),'m');

