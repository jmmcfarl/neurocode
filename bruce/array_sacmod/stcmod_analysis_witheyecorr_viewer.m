clear all
close all

dir_prefix = '~';
Expt_nums = [86 93];

all_ov_info = [];
all_gsac_info = [];
all_msac_info = [];
all_avgrate = [];
all_nspks = [];
all_gsac_rate = [];
all_msac_rate = [];
for ex = 1:length(Expt_nums)
    
    
    load ~/Analysis/bruce/summary_analysis/su_data.mat
    Expt_name = sprintf('G0%d',Expt_nums(ex));
    data_dir = [dir_prefix '/Data/bruce/' Expt_name];
    anal_dir = ['~/Analysis/bruce/' Expt_name];
    
    %%
    cd(anal_dir)
    anal_name = 'hbar_stc_sacmod_analysis';
    load(anal_name);
    
    sua_avgrate_hbar = [sua_data(:).avg_rate];
    sua_nspkps_hbar = [sua_data(:).nspks];
    
    sua_gsac_info_hbar = [sua_data(:).gsacdep_info];
    sua_ovg_info_hbar = [sua_data(:).ov_info_gsac];
    sua_msac_info_hbar = [sua_data(:).msacdep_info];
    sua_ovm_info_hbar = [sua_data(:).ov_info_msac];
    
    sac_tick = sua_data(1).Stick;
    sua_gmarg_sacrate_hbar = reshape([sua_data(:).marg_gsacrate],length(sac_tick),length(sua_data));
    sua_mmarg_sacrate_hbar = reshape([sua_data(:).marg_msacrate],length(sac_tick),length(sua_data));
    
    sua_data_h = sua_data;
    mua_data_h = mua_data;
    %%
    cd(anal_dir)
    anal_name = 'vbar_stc_sacmod_analysis';
    load(anal_name);
    
    sua_avgrate_vbar = [sua_data(:).avg_rate];
    sua_nspkps_vbar = [sua_data(:).nspks];
    
    sua_gsac_info_vbar = [sua_data(:).gsacdep_info];
    sua_ovg_info_vbar = [sua_data(:).ov_info_gsac];
    sua_msac_info_vbar = [sua_data(:).msacdep_info];
    sua_ovm_info_vbar = [sua_data(:).ov_info_msac];
    
    sua_gmarg_sacrate_vbar = reshape([sua_data(:).marg_gsacrate],length(sac_tick),length(sua_data));
    sua_mmarg_sacrate_vbar = reshape([sua_data(:).marg_msacrate],length(sac_tick),length(sua_data));
    
    sua_data_v = sua_data;
    mua_data_v = mua_data;
    %%
    all_ov_info = cat(1,all_ov_info,[sua_ovg_info_hbar(:) sua_ovg_info_vbar(:)]);
    all_avgrate = cat(1,all_avgrate,[sua_avgrate_hbar(:) sua_avgrate_vbar(:)]);
    all_nspks = cat(1,all_nspks,[sua_nspkps_hbar(:) sua_nspkps_vbar(:)]);
    
    temp = cat(3,sua_gsac_info_hbar,sua_gsac_info_vbar);
    all_gsac_info = cat(2,all_gsac_info,temp);
    temp = cat(3,sua_gmarg_sacrate_hbar,sua_gmarg_sacrate_vbar);
    all_gsac_rate = cat(2,all_gsac_rate,temp);
    
    temp = cat(3,sua_msac_info_hbar,sua_msac_info_vbar);
    all_msac_info = cat(2,all_msac_info,temp);
    temp = cat(3,sua_mmarg_sacrate_hbar,sua_mmarg_sacrate_vbar);
    all_msac_rate = cat(2,all_msac_rate,temp);
end

%% FOR A SINGLE MUA EXPT
for ss = 1:96
   if ~isempty(mua_data_h(ss).avg_rate)
       ov_info_h(ss) = mua_data_h(ss).ov_info_gsac;
       norm_grate_h(ss,:) = mua_data_h(ss).marg_gsacrate/mua_data_h(ss).avg_rate;
       norm_ginfo_h(ss,:) = mua_data_h(ss).gsacdep_info/mua_data_h(ss).ov_info_gsac;
       ginfo_spk = mua_data_h(ss).gsacdep_info./mua_data_h(ss).marg_gsacrate';
       ov_info_spk_h(ss) = ov_info_h(ss)/mua_data_h(ss).avg_rate;
       norm_ginfo_spk_h(ss,:) = ginfo_spk/ov_info_spk_h(ss);

       norm_mrate_h(ss,:) = mua_data_h(ss).marg_msacrate/mua_data_h(ss).avg_rate;
       norm_minfo_h(ss,:) = mua_data_h(ss).msacdep_info/mua_data_h(ss).ov_info_msac;
       minfo_spk = mua_data_h(ss).msacdep_info./mua_data_h(ss).marg_msacrate';
       norm_minfo_spk_h(ss,:) = minfo_spk/ov_info_spk_h(ss);
   else
       ov_info_h(ss) = nan;
       ov_info_spk_h(ss) = nan;
       norm_grate_h(ss,:) = nan;
       norm_ginfo_h(ss,:) = nan;
       norm_ginfo_spk_h(ss,:) = nan;
       norm_mrate_h(ss,:) = nan;
       norm_minfo_h(ss,:) = nan;
       norm_minfo_spk_h(ss,:) = nan;
   end
   
   if ~isempty(mua_data_v(ss).avg_rate)
       ov_info_v(ss) = mua_data_v(ss).ov_info_gsac;
       
       norm_grate_v(ss,:) = mua_data_v(ss).marg_gsacrate/mua_data_v(ss).avg_rate;
       norm_ginfo_v(ss,:) = mua_data_v(ss).gsacdep_info/mua_data_v(ss).ov_info_gsac;
       ginfo_spk = mua_data_v(ss).gsacdep_info./mua_data_v(ss).marg_gsacrate';
       ov_info_spk_v(ss) = ov_info_v(ss)/mua_data_v(ss).avg_rate;
       norm_ginfo_spk_v(ss,:) = ginfo_spk/ov_info_spk_v(ss);
       
       norm_mrate_v(ss,:) = mua_data_v(ss).marg_msacrate/mua_data_v(ss).avg_rate;
       norm_minfo_v(ss,:) = mua_data_v(ss).msacdep_info/mua_data_v(ss).ov_info_msac;
       minfo_spk = mua_data_v(ss).msacdep_info./mua_data_v(ss).marg_msacrate';
       norm_minfo_spk_v(ss,:) = minfo_spk/ov_info_spk_v(ss);
   else
       ov_info_v(ss) = nan;
       ov_info_spk_v(ss) = nan;
       norm_grate_v(ss,:) = nan;
       norm_ginfo_v(ss,:) = nan;
       norm_ginfo_spk_v(ss,:) = nan;
       norm_mrate_v(ss,:) = nan;
       norm_minfo_v(ss,:) = nan;
       norm_minfo_spk_v(ss,:) = nan;
   end    
   
   if ~isempty(mua_data_v(ss).avg_rate) && ~isempty(mua_data_h(ss).avg_rate)
       if ov_info_h(ss) > ov_info_v(ss)
           ov_info_b(ss) = ov_info_h(ss);
           ov_info_spk_b(ss) = ov_info_spk_h(ss);
           norm_grate_b(ss,:) = norm_grate_h(ss,:);
           norm_ginfo_b(ss,:) = norm_ginfo_h(ss,:);
           norm_ginfo_spk_b(ss,:) = norm_ginfo_spk_h(ss,:);
            norm_mrate_b(ss,:) = norm_mrate_h(ss,:);
           norm_minfo_b(ss,:) = norm_minfo_h(ss,:);
           norm_minfo_spk_b(ss,:) = norm_minfo_spk_h(ss,:);          
       else
           ov_info_b(ss) = ov_info_v(ss);
           ov_info_spk_b(ss) = ov_info_spk_v(ss);
           norm_grate_b(ss,:) = norm_grate_v(ss,:);
           norm_ginfo_b(ss,:) = norm_ginfo_v(ss,:);
           norm_ginfo_spk_b(ss,:) = norm_ginfo_spk_v(ss,:);
            norm_mrate_b(ss,:) = norm_mrate_v(ss,:);
           norm_minfo_b(ss,:) = norm_minfo_v(ss,:);
           norm_minfo_spk_b(ss,:) = norm_minfo_spk_v(ss,:);                     
       end
   end
end

Stick = mua_data(1).Stick;
xb = 2;xe = length(Stick)-1;

close all
min_info = 0.02;
use_b = find(ov_info_b > min_info);
figure;
h1=shadedErrorBar(Stick,mean(norm_grate_b(use_h,:)),std(norm_grate_b(use_h,:))/sqrt(length(use_b)),{'color','b'});
hold on
h2=shadedErrorBar(Stick,mean(norm_mrate_b(use_v,:)),std(norm_mrate_b(use_v,:))/sqrt(length(use_b)),{'color','r'});
% h1=shadedErrorBar(Stick,mean(norm_grate_b(use_h,:)),std(norm_grate_b(use_h,:)));
% hold on
% h2=shadedErrorBar(Stick,mean(norm_mrate_b(use_v,:)),std(norm_mrate_b(use_v,:)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
xlim(Stick([xb xe]));
xlabel('Time since saccade onset (s)','fontsize',12)
ylabel('Relative firing rate','fontsize',12);
xl = xlim();
line(xl,[1 1],'color','k');
box off;
fillPage(gcf,'papersize',[5 5]);

figure;
h1=shadedErrorBar(Stick,mean(norm_ginfo_b(use_h,:)),std(norm_ginfo_b(use_h,:))/sqrt(length(use_b)),{'color','b'});
hold on
h2=shadedErrorBar(Stick,mean(norm_minfo_b(use_v,:)),std(norm_minfo_b(use_v,:))/sqrt(length(use_b)),{'color','r'});
% h1=shadedErrorBar(Stick,mean(norm_ginfo_b(use_h,:)),std(norm_ginfo_b(use_h,:)));
% hold on
% h2=shadedErrorBar(Stick,mean(norm_minfo_b(use_v,:)),std(norm_minfo_b(use_v,:)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
xlim(Stick([xb xe]));
xlabel('Time since saccade onset (s)','fontsize',12)
ylabel('Relative info (per sample)','fontsize',12);
xl = xlim();
line(xl,[1 1],'color','k');
box off;
fillPage(gcf,'papersize',[5 5]);

figure;
h1=shadedErrorBar(Stick,mean(norm_ginfo_spk_b(use_h,:)),std(norm_ginfo_spk_b(use_h,:))/sqrt(length(use_b)),{'color','b'});
hold on
h2=shadedErrorBar(Stick,mean(norm_minfo_spk_b(use_v,:)),std(norm_minfo_spk_b(use_v,:))/sqrt(length(use_b)),{'color','r'});
% h1=shadedErrorBar(Stick,mean(norm_ginfo_spk_b(use_h,:)),std(norm_ginfo_spk_b(use_h,:)));
% hold on
% h2=shadedErrorBar(Stick,mean(norm_minfo_spk_b(use_v,:)),std(norm_minfo_spk_b(use_v,:)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
xlim(Stick([xb xe]));
xlabel('Time since saccade onset (s)','fontsize',12)
ylabel('Relative info (per spike)','fontsize',12);
xl = xlim();
line(xl,[1 1],'color','k');
box off;
fillPage(gcf,'papersize',[5 5]);

% use_h = find(~isnan(ov_info_h));
% use_v = find(~isnan(ov_info_v));
% use_h = use_h(ov_info_h(use_h) > min_info);
% use_v = use_v(ov_info_v(use_v) > min_info);
% figure;
% shadedErrorBar(Stick,mean(norm_grate_h(use_h,:)),std(norm_grate_h(use_h,:))/sqrt(length(use_h)));
% hold on
% shadedErrorBar(Stick,mean(norm_grate_v(use_v,:)),std(norm_grate_v(use_v,:))/sqrt(length(use_v)),{'color','r'});
% 
% figure;
% shadedErrorBar(Stick,mean(norm_ginfo_h(use_h,:)),std(norm_ginfo_h(use_h,:))/sqrt(length(use_h)));
% hold on
% shadedErrorBar(Stick,mean(norm_ginfo_v(use_v,:)),std(norm_ginfo_v(use_v,:))/sqrt(length(use_v)),{'color','r'});
% 
% figure;
% shadedErrorBar(Stick,mean(norm_ginfo_spk_h(use_h,:)),std(norm_ginfo_spk_h(use_h,:))/sqrt(length(use_h)));
% hold on
% shadedErrorBar(Stick,mean(norm_ginfo_spk_v(use_v,:)),std(norm_ginfo_spk_v(use_v,:))/sqrt(length(use_v)),{'color','r'});

%%
min_nspks = 5e3;
min_tune = 1.5;

norm_grate = bsxfun(@rdivide,all_gsac_rate,reshape(all_avgrate,[1 size(all_avgrate)]));
% [best_rate,pref_ori] = max(all_avgrate,[],2);
[best_info,pref_ori] = max(all_ov_info,[],2);
norm_ginfo = bsxfun(@rdivide,all_gsac_info,reshape(all_ov_info,[1 size(all_ov_info)]));
norm_mrate = bsxfun(@rdivide,all_msac_rate,reshape(all_avgrate,[1 size(all_avgrate)]));
norm_minfo = bsxfun(@rdivide,all_msac_info,reshape(all_ov_info,[1 size(all_avgrate)]));
for ii = 1:size(all_avgrate,1)
    nonpref_ori(ii) = setdiff([1 2],pref_ori(ii));
    pref_nspks(ii) = all_nspks(ii,pref_ori(ii));
    npref_nspks(ii) = all_nspks(ii,nonpref_ori(ii));
    pref_rate(ii) = all_avgrate(ii,pref_ori(ii));
    npref_rate(ii) = all_avgrate(ii,nonpref_ori(ii));
    pref_grate(ii,:) = squeeze(norm_grate(:,ii,pref_ori(ii)));
    npref_grate(ii,:) = squeeze(norm_grate(:,ii,nonpref_ori(ii)));
    pref_ginfo(ii,:) = squeeze(norm_ginfo(:,ii,pref_ori(ii)));
    npref_ginfo(ii,:) = squeeze(norm_ginfo(:,ii,nonpref_ori(ii)));
    pref_ovinfo(ii) = all_ov_info(ii,pref_ori(ii));
    npref_ovinfo(ii) = all_ov_info(ii,nonpref_ori(ii));
    pref_mrate(ii,:) = squeeze(norm_mrate(:,ii,pref_ori(ii)));
    npref_mrate(ii,:) = squeeze(norm_mrate(:,ii,nonpref_ori(ii)));
    pref_minfo(ii,:) = squeeze(norm_minfo(:,ii,pref_ori(ii)));
    npref_minfo(ii,:) = squeeze(norm_minfo(:,ii,nonpref_ori(ii)));
end
% usable = find(min(all_nspks,[],2) > min_nspks);
% ori_tune = pref_ovinfo./npref_ovinfo;
% usable = usable(ori_tune(usable) > min_tune);

% close all
% Stick = sua_data(1).Stick;
% figure;
% shadedErrorBar(Stick,mean(pref_grate(usable,:)),std(pref_grate(usable,:))/sqrt(length(usable)));
% hold on
% shadedErrorBar(Stick,mean(npref_grate(usable,:)),std(npref_grate(usable,:))/sqrt(length(usable)),{'color','r'});
%
% figure;
% shadedErrorBar(Stick,mean(pref_ginfo(usable,:)),std(pref_ginfo(usable,:))/sqrt(length(usable)));
% hold on
% shadedErrorBar(Stick,mean(npref_ginfo(usable,:)),std(npref_ginfo(usable,:))/sqrt(length(usable)),{'color','r'});

close all
% min_info = 0.2;
% info_perspk = best_info./pref_rate';
% usable = find(pref_nspks > min_nspks & info_perspk' > min_info);
min_info = 0.02;
usable = find(pref_nspks > min_nspks & best_info' > min_info);
figure;
h1=shadedErrorBar(Stick,mean(pref_grate(usable,:)),std(pref_grate(usable,:))/sqrt(length(usable)),{'color','b'});
hold on
h2=shadedErrorBar(Stick,mean(pref_mrate(usable,:)),std(pref_mrate(usable,:))/sqrt(length(usable)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
xlabel('Time since saccade onset (s)','fontsize',12)
ylabel('Average relative rate','fontsize',12);
xlim([-0.2 0.5]);
xl = xlim();
line(xl,[1 1],'color','k');
box off;
fillPage(gcf,'papersize',[5 5]);

figure;
h1=shadedErrorBar(Stick,mean(pref_ginfo(usable,:)),std(pref_ginfo(usable,:))/sqrt(length(usable)),{'color','b'});
hold on
h2=shadedErrorBar(Stick,mean(pref_minfo(usable,:)),std(pref_minfo(usable,:))/sqrt(length(usable)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
xlabel('Time since saccade onset (s)','fontsize',12)
ylabel('Relative information (per sample)','fontsize',12);
xlim([-0.2 0.5]);
xl = xlim();
line(xl,[1 1],'color','k');
box off;
fillPage(gcf,'papersize',[5 5]);

figure;
h1=shadedErrorBar(Stick,mean(pref_ginfo(usable,:)./pref_grate(usable,:)),std(pref_ginfo(usable,:)./pref_grate(usable,:))/sqrt(length(usable)),{'color','b'});
hold on
h2=shadedErrorBar(Stick,mean(pref_minfo(usable,:)./pref_mrate(usable,:)),std(pref_minfo(usable,:)./pref_mrate(usable,:))/sqrt(length(usable)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
xlabel('Time since saccade onset (s)','fontsize',12)
ylabel('Relative information (per spike)','fontsize',12);
xlim([-0.2 0.5]);
xl = xlim();
line(xl,[1 1],'color','k');
box off;
fillPage(gcf,'papersize',[5 5]);

% figure;
% h1=shadedErrorBar(Stick,mean(npref_grate(usable,:)),std(npref_grate(usable,:))/sqrt(length(usable)));
% hold on
% h2=shadedErrorBar(Stick,mean(npref_mrate(usable,:)),std(npref_mrate(usable,:))/sqrt(length(usable)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
% xlabel('Time since saccade onset (s)','fontsize',12)
% ylabel('Average relative rate','fontsize',12);
% xlim([-0.2 0.5]);
% xl = xlim();
% line(xl,[1 1],'color','k');
%
% figure;
% h1=shadedErrorBar(Stick,mean(npref_ginfo(usable,:)),std(npref_ginfo(usable,:))/sqrt(length(usable)));
% hold on
% h2=shadedErrorBar(Stick,mean(npref_minfo(usable,:)),std(npref_minfo(usable,:))/sqrt(length(usable)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
% xlabel('Time since saccade onset (s)','fontsize',12)
% ylabel('Relative information (per sample)','fontsize',12);
% xlim([-0.2 0.5]);
% xl = xlim();
% line(xl,[1 1],'color','k');
%
% figure;
% shadedErrorBar(Stick,mean(npref_ginfo(usable,:)./npref_grate(usable,:)),std(npref_ginfo(usable,:)./npref_grate(usable,:))/sqrt(length(usable)));
% hold on
% shadedErrorBar(Stick,mean(npref_minfo(usable,:)./npref_mrate(usable,:)),std(npref_minfo(usable,:)./npref_mrate(usable,:))/sqrt(length(usable)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'Guided sacs','Micro sacs'});
% xlabel('Time since saccade onset (s)','fontsize',12)
% ylabel('Relative information (per spike)','fontsize',12);
% xlim([-0.2 0.5]);
% xl = xlim();
% line(xl,[1 1],'color','k');

%% FOR SUA
flen = 12;
use_nPix = 24;
dt = 0.01;
Xtick = sua_data(1).Stick;
xl = sac_tick([1 end]);
close all
f1 = figure('Name','Horizontal');
f2 = figure('Name','Vertical');
f3 = figure('Name','Hori mods');
f4 = figure('Name','Vert mods');
for ss = 1:length(sua_data)
    
    fprintf('Viewing SU %d of %d\n',ss,length(sua_data_h));
    fprintf('Hor: rate: %.2f  nspks: %d\n',sua_avgrate_hbar(ss)/dt,sua_nspkps_hbar(ss));
    fprintf('Ver: rate: %.2f  nspks: %d\n',sua_avgrate_vbar(ss)/dt,sua_nspkps_vbar(ss));
    
    marg_gdist = sua_data_h(ss).marg_gdist_gsac;
    cum_dist = cumsum(marg_gdist);
    yb = find(cum_dist > 0.001,1,'first');
    ye = find(cum_dist > 0.995,1,'first');
    xb = 2;xe = length(Xtick)-1;
    Ytick = sua_data_h(ss).Gtick;
    marg_gsacrate = sua_data_h(ss).marg_gsacrate/dt;
    zp = find(Xtick >= 0,1);
    ze = find(Xtick >= 0.2,1);
    [~,temp] = min(marg_gsacrate(zp:ze));
    mloc = zp + temp-1;
    [~,temp] = max(marg_gsacrate(zp:ze));
    mxloc = zp + temp - 1;
    
    marg_grate = sum(sua_data_h(ss).gsac_TB_dist.*sua_data_h(ss).gsac_TB_rate,2)./marg_gdist;
    
    % PLOT GSAC
    figure(f1);clf
    subplot(4,3,[1 2])
    pcolor(Xtick,Ytick(yb:ye),sua_data_h(ss).gsac_TB_rate(yb:ye,:)/dt);shading flat
    ca = caxis(); caxis([ca(1) ca(end)*0.8])
    % colorbar
    ylim(Ytick([yb ye]));
    xlim(Xtick([xb xe]));
    line(Xtick([mloc mloc]),Ytick([yb ye]),'color','r');
    line(Xtick([mxloc mxloc]),Ytick([yb ye]),'color','g');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Generating signal','fontsize',14)
    subplot(4,3,[3])
    plot(marg_gdist,Ytick,'b'); set(gca,'xtick',[],'ytick',[]);
    ylim(Ytick([yb ye]));
    
    subplot(4,3,[4 5])
    plot(Xtick,marg_gsacrate);hold on
    yl = ylim();
    line(Xtick([mloc mloc]),yl,'color','r');
    line(Xtick([mxloc mxloc]),yl,'color','g');
    xlim(Xtick([xb xe]));
    line(Xtick([xb xe]),ones(1,2)*sua_data_h(ss).avg_rate/dt,'color','k');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Firing rate (Hz)','fontsize',14)
    
    subplot(4,3,6)
    hold on
    plot(Ytick,sua_data_h(ss).gsac_TB_rate(:,mloc)/dt,'r')
    plot(Ytick,sua_data_h(ss).gsac_TB_rate(:,mxloc)/dt,'g')
    plot(Ytick,marg_grate/dt,'k')
    yl = ylim();
    plot(Ytick,marg_gdist/max(marg_gdist)*range(yl)*0.8,'b')
    xlim(Ytick([yb ye]));
    xlabel('Generating signal','fontsize',14)
    ylabel('Firing rate (Hz)','fontsize',14)
    
    subplot(4,3,[7 8])
    plot(Xtick,sua_data_h(ss).gsacdep_info/dt);hold on
    xlim(Xtick([xb xe]));
    line(Xtick([xb xe]),ones(1,2)*sua_data_h(ss).ov_info_gsac/dt,'color','k')
    yl = ylim(); ylim([0 yl(2)]); yl = ylim();
    line(Xtick([mloc mloc]),yl,'color','r');
    line(Xtick([mxloc mxloc]),yl,'color','g');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Stimulus info (bits/sec)','fontsize',14)
    
    subplot(4,3,[10 11])
    plot(Xtick,sua_data_h(ss).gsacdep_info./marg_gsacrate'/dt);hold on
    xlim(Xtick([xb xe]));
    line(Xtick([xb xe]),ones(1,2)*sua_data_h(ss).ov_info_gsac/sua_data_h(ss).avg_rate,'color','k')
    yl = ylim(); ylim([0 yl(2)]); yl = ylim();
    line(Xtick([mloc mloc]),yl,'color','r');
    line(Xtick([mxloc mxloc]),yl,'color','g');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Stimulus info (bits/spk)','fontsize',14)
    
    %
    marg_gdist = sua_data_v(ss).marg_gdist_gsac;
    cum_dist = cumsum(marg_gdist);
    yb = find(cum_dist > 0.001,1,'first');
    ye = find(cum_dist > 0.995,1,'first');
    xb = 2;xe = length(Xtick)-1;
    Ytick = sua_data_v(ss).Gtick;
    marg_gsacrate = sua_data_v(ss).marg_gsacrate'/dt;
    zp = find(Xtick >= 0,1);
    ze = find(Xtick >= 0.2,1);
    [~,temp] = min(marg_gsacrate(zp:ze));
    mloc = zp + temp-1;
    [~,temp] = max(marg_gsacrate(zp:ze));
    mxloc = zp + temp - 1;
    
    marg_grate = sum(sua_data_v(ss).gsac_TB_dist.*sua_data_v(ss).gsac_TB_rate,2)./marg_gdist;
    
    % PLOT GSAC
    figure(f2);clf
    subplot(4,3,[1 2])
    pcolor(Xtick,Ytick(yb:ye),sua_data_v(ss).gsac_TB_rate(yb:ye,:)/dt);shading flat
    % colorbar
    ca = caxis(); caxis([ca(1) ca(end)*0.8])
    ylim(Ytick([yb ye]));
    xlim(Xtick([xb xe]));
    line(Xtick([mloc mloc]),Ytick([yb ye]),'color','r');
    line(Xtick([mxloc mxloc]),Ytick([yb ye]),'color','g');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Generating signal','fontsize',14)
    subplot(4,3,[3])
    plot(marg_gdist,Ytick,'b'); set(gca,'xtick',[],'ytick',[]);
    ylim(Ytick([yb ye]));
    
    subplot(4,3,[4 5])
    plot(Xtick,marg_gsacrate);hold on
    yl = ylim();
    line(Xtick([mloc mloc]),yl,'color','r');
    line(Xtick([mxloc mxloc]),yl,'color','g');
    xlim(Xtick([xb xe]));
    line(Xtick([xb xe]),ones(1,2)*sua_data_v(ss).avg_rate/dt,'color','k');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Firing rate (Hz)','fontsize',14)
    
    subplot(4,3,6)
    hold on
    plot(Ytick,sua_data_v(ss).gsac_TB_rate(:,mloc)/dt,'r')
    plot(Ytick,sua_data_v(ss).gsac_TB_rate(:,mxloc)/dt,'g')
    plot(Ytick,marg_grate/dt,'k')
    yl = ylim();
    plot(Ytick,marg_gdist/max(marg_gdist)*range(yl)*0.8,'b')
    xlim(Ytick([yb ye]));
    xlabel('Generating signal','fontsize',14)
    ylabel('Firing rate (Hz)','fontsize',14)
    
    subplot(4,3,[7 8])
    plot(Xtick,sua_data_v(ss).gsacdep_info/dt);hold on
    xlim(Xtick([xb xe]));
    line(Xtick([xb xe]),ones(1,2)*sua_data_v(ss).ov_info_gsac/dt,'color','k')
    yl = ylim(); ylim([0 yl(2)]); yl = ylim();
    line(Xtick([mloc mloc]),yl,'color','r');
    line(Xtick([mxloc mxloc]),yl,'color','g');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Stimulus info (bits/sec)','fontsize',14)
    
    subplot(4,3,[10 11])
    plot(Xtick,sua_data_v(ss).gsacdep_info./marg_gsacrate/dt);hold on
    xlim(Xtick([xb xe]));
    line(Xtick([xb xe]),ones(1,2)*sua_data_v(ss).ov_info_gsac/sua_data_v(ss).avg_rate,'color','k')
    yl = ylim(); ylim([0 yl(2)]); yl = ylim();
    line(Xtick([mloc mloc]),yl,'color','r');
    line(Xtick([mxloc mxloc]),yl,'color','g');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Stimulus info (bits/spk)','fontsize',14)
    
    
    %
    figure(f3); clf
    cur_sta = sua_data_h(ss).sta; ca = max(abs(cur_sta));
    subplot(3,3,1)
    imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
    title('STA')
    cur_stcs = sua_data_h(ss).stcs; ca = max(abs(cur_stcs(:)));
    for ii = 1:3
        subplot(3,3,3+ii)
        imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
        title(sprintf('Exc STC %d',ii));
    end
    for ii = 1:3
        subplot(3,3,6+ii)
        imagesc(reshape(cur_stcs(:,end-3+ii),flen,use_nPix)); caxis([-ca ca]);
         title(sprintf('Sup STC %d',ii));
   end
    
    %
    figure(f4); clf
    cur_sta = sua_data_v(ss).sta; ca = max(abs(cur_sta));
    subplot(3,3,1)
    imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
    title('STA')
    cur_stcs = sua_data_v(ss).stcs; ca = max(abs(cur_stcs(:)));
    for ii = 1:3
        subplot(3,3,3+ii)
        imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
        title(sprintf('Exc STC %d',ii));
    end
    for ii = 1:3
        subplot(3,3,6+ii)
        imagesc(reshape(cur_stcs(:,end-3+ii),flen,use_nPix)); caxis([-ca ca]);
         title(sprintf('Sup STC %d',ii));
    end
    
    pause
end


%% FOR MUA
flen = 12;
use_nPix = 24;
dt = 0.01;
Xtick = mua_data(1).Stick;
xl = sac_tick([1 end]);
close all
f1 = figure('Name','Horizontal');
f2 = figure('Name','Vertical');
f3 = figure('Name','Hori mods');
f4 = figure('Name','Vert mods');
for ss = 1:length(mua_data)
    
    if ~isempty(mua_data_h(ss).nspks)
        fprintf('Viewing SU %d of %d\n',ss,length(mua_data_h));
        fprintf('Hor: rate: %.2f  nspks: %d\n',mua_data_h(ss).avg_rate/dt,mua_data_h(ss).nspks);
        fprintf('Ver: rate: %.2f  nspks: %d\n',mua_data_v(ss).avg_rate/dt,mua_data_v(ss).nspks);
        
        marg_gdist = mua_data_h(ss).marg_gdist_gsac;
        cum_dist = cumsum(marg_gdist);
        yb = find(cum_dist > 0.001,1,'first');
        ye = find(cum_dist > 0.995,1,'first');
        xb = 2;xe = length(Xtick)-1;
        Ytick = mua_data_h(ss).Gtick;
        marg_gsacrate = mua_data_h(ss).marg_gsacrate/dt;
        zp = find(Xtick >= 0,1);
        ze = find(Xtick >= 0.2,1);
        [~,temp] = min(marg_gsacrate(zp:ze));
        mloc = zp + temp-1;
        [~,temp] = max(marg_gsacrate(zp:ze));
        mxloc = zp + temp - 1;
        
        marg_grate = sum(mua_data_h(ss).gsac_TB_dist.*mua_data_h(ss).gsac_TB_rate,2)./marg_gdist;
        
        % PLOT GSAC
        figure(f1);clf
        subplot(4,3,[1 2])
        pcolor(Xtick,Ytick(yb:ye),mua_data_h(ss).gsac_TB_rate(yb:ye,:)/dt);shading flat
        % colorbar
        ylim(Ytick([yb ye]));
        xlim(Xtick([xb xe]));
        line(Xtick([mloc mloc]),Ytick([yb ye]),'color','r');
        line(Xtick([mxloc mxloc]),Ytick([yb ye]),'color','g');
        xlabel('Time since saccade onset (s)','fontsize',14)
        ylabel('Generating signal','fontsize',14)
        subplot(4,3,[3])
        plot(marg_gdist,Ytick,'b'); set(gca,'xtick',[],'ytick',[]);
        ylim(Ytick([yb ye]));
        
        subplot(4,3,[4 5])
        plot(Xtick,marg_gsacrate);hold on
        yl = ylim();
        line(Xtick([mloc mloc]),yl,'color','r');
        line(Xtick([mxloc mxloc]),yl,'color','g');
        xlim(Xtick([xb xe]));
        line(Xtick([xb xe]),ones(1,2)*mua_data_h(ss).avg_rate/dt,'color','k');
        xlabel('Time since saccade onset (s)','fontsize',14)
        ylabel('Firing rate (Hz)','fontsize',14)
        
        subplot(4,3,6)
        hold on
        plot(Ytick,mua_data_h(ss).gsac_TB_rate(:,mloc)/dt,'r')
        plot(Ytick,mua_data_h(ss).gsac_TB_rate(:,mxloc)/dt,'g')
        plot(Ytick,marg_grate/dt,'k')
        yl = ylim();
        plot(Ytick,marg_gdist/max(marg_gdist)*range(yl)*0.8,'b')
        xlim(Ytick([yb ye]));
        xlabel('Generating signal','fontsize',14)
        ylabel('Firing rate (Hz)','fontsize',14)
        
        subplot(4,3,[7 8])
        plot(Xtick,mua_data_h(ss).gsacdep_info/dt);hold on
        xlim(Xtick([xb xe]));
        line(Xtick([xb xe]),ones(1,2)*mua_data_h(ss).ov_info_gsac/dt,'color','k')
        yl = ylim(); ylim([0 yl(2)]); yl = ylim();
        line(Xtick([mloc mloc]),yl,'color','r');
        line(Xtick([mxloc mxloc]),yl,'color','g');
        xlabel('Time since saccade onset (s)','fontsize',14)
        ylabel('Stimulus info (bits/sec)','fontsize',14)
        
        subplot(4,3,[10 11])
        plot(Xtick,mua_data_h(ss).gsacdep_info./marg_gsacrate'/dt);hold on
        xlim(Xtick([xb xe]));
        line(Xtick([xb xe]),ones(1,2)*mua_data_h(ss).ov_info_gsac/mua_data_h(ss).avg_rate,'color','k')
        yl = ylim(); ylim([0 yl(2)]); yl = ylim();
        line(Xtick([mloc mloc]),yl,'color','r');
        line(Xtick([mxloc mxloc]),yl,'color','g');
        xlabel('Time since saccade onset (s)','fontsize',14)
        ylabel('Stimulus info (bits/spk)','fontsize',14)
        
        
        figure(f3); clf
        cur_sta = mua_data_h(ss).sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
        cur_stcs = mua_data_h(ss).stcs; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
        end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-3+ii),flen,use_nPix)); caxis([-ca ca]);
        end
        
    end
    %%
    if ~isempty(mua_data_v(ss).nspks)
        
        marg_gdist = mua_data_v(ss).marg_gdist_gsac;
        cum_dist = cumsum(marg_gdist);
        yb = find(cum_dist > 0.001,1,'first');
        ye = find(cum_dist > 0.995,1,'first');
        xb = 2;xe = length(Xtick)-1;
        Ytick = mua_data_v(ss).Gtick;
        marg_gsacrate = mua_data_v(ss).marg_gsacrate'/dt;
        zp = find(Xtick >= 0,1);
        ze = find(Xtick >= 0.2,1);
        [~,temp] = min(marg_gsacrate(zp:ze));
        mloc = zp + temp-1;
        [~,temp] = max(marg_gsacrate(zp:ze));
        mxloc = zp + temp - 1;
        
        marg_grate = sum(mua_data_v(ss).gsac_TB_dist.*mua_data_v(ss).gsac_TB_rate,2)./marg_gdist;
        
        % PLOT GSAC
        figure(f2);clf
        subplot(4,3,[1 2])
        pcolor(Xtick,Ytick(yb:ye),mua_data_v(ss).gsac_TB_rate(yb:ye,:)/dt);shading flat
        % colorbar
        ylim(Ytick([yb ye]));
        xlim(Xtick([xb xe]));
        line(Xtick([mloc mloc]),Ytick([yb ye]),'color','r');
        line(Xtick([mxloc mxloc]),Ytick([yb ye]),'color','g');
        xlabel('Time since saccade onset (s)','fontsize',14)
        ylabel('Generating signal','fontsize',14)
        subplot(4,3,[3])
        plot(marg_gdist,Ytick,'b'); set(gca,'xtick',[],'ytick',[]);
        ylim(Ytick([yb ye]));
        
        subplot(4,3,[4 5])
        plot(Xtick,marg_gsacrate);hold on
        yl = ylim();
        line(Xtick([mloc mloc]),yl,'color','r');
        line(Xtick([mxloc mxloc]),yl,'color','g');
        xlim(Xtick([xb xe]));
        line(Xtick([xb xe]),ones(1,2)*mua_data_v(ss).avg_rate/dt,'color','k');
        xlabel('Time since saccade onset (s)','fontsize',14)
        ylabel('Firing rate (Hz)','fontsize',14)
        
        subplot(4,3,6)
        hold on
        plot(Ytick,mua_data_v(ss).gsac_TB_rate(:,mloc)/dt,'r')
        plot(Ytick,mua_data_v(ss).gsac_TB_rate(:,mxloc)/dt,'g')
        plot(Ytick,marg_grate/dt,'k')
        yl = ylim();
        plot(Ytick,marg_gdist/max(marg_gdist)*range(yl)*0.8,'b')
        xlim(Ytick([yb ye]));
        xlabel('Generating signal','fontsize',14)
        ylabel('Firing rate (Hz)','fontsize',14)
        
        subplot(4,3,[7 8])
        plot(Xtick,mua_data_v(ss).gsacdep_info/dt);hold on
        xlim(Xtick([xb xe]));
        line(Xtick([xb xe]),ones(1,2)*mua_data_v(ss).ov_info_gsac/dt,'color','k')
        yl = ylim(); ylim([0 yl(2)]); yl = ylim();
        line(Xtick([mloc mloc]),yl,'color','r');
        line(Xtick([mxloc mxloc]),yl,'color','g');
        xlabel('Time since saccade onset (s)','fontsize',14)
        ylabel('Stimulus info (bits/sec)','fontsize',14)
        
        subplot(4,3,[10 11])
        plot(Xtick,mua_data_v(ss).gsacdep_info./marg_gsacrate/dt);hold on
        xlim(Xtick([xb xe]));
        line(Xtick([xb xe]),ones(1,2)*mua_data_v(ss).ov_info_gsac/mua_data_v(ss).avg_rate,'color','k')
        yl = ylim(); ylim([0 yl(2)]); yl = ylim();
        line(Xtick([mloc mloc]),yl,'color','r');
        line(Xtick([mxloc mxloc]),yl,'color','g');
        xlabel('Time since saccade onset (s)','fontsize',14)
        ylabel('Stimulus info (bits/spk)','fontsize',14)
        
        %%
        
        %%
        figure(f4); clf
        cur_sta = mua_data_v(ss).sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
        cur_stcs = mua_data_v(ss).stcs; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
        end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-3+ii),flen,use_nPix)); caxis([-ca ca]);
        end
    end
    
    pause
end