clear all
close all


poss_bar_oris = [0 45 90 135];
for EE = 1:4
    sname = sprintf('unit_sac_kern_mods_unsm_%ddeg',poss_bar_oris(EE));
    load(sname);
    
    gray_sac_avg(EE,:,:) = gray_sac_kern;
    gray_sac_sem(EE,:,:) = gray_sac_se;
    gray_osac_avg(EE,:,:) = gray_osac_kern2;
    gray_osac_sem(EE,:,:) = gray_osac_se2;
    gray_rsac_avg(EE,:,:) = gray_rsac_kern2;
    gray_rsac_sem(EE,:,:) = gray_rsac_se2;
    gray_msac_avg(EE,:,:) = gray_msac_kern2;
    gray_msac_sem(EE,:,:) = gray_msac_se2;
 
    im_sac_avg(EE,:,:) = im_sac_kern;
    im_sac_sem(EE,:,:) = im_sac_se;
    im_osac_avg(EE,:,:) = im_rsac_kern2; %o and r sacs were switched 
    im_osac_sem(EE,:,:) = im_rsac_se2;
    im_rsac_avg(EE,:,:) = im_osac_kern2;
    im_rsac_sem(EE,:,:) = im_osac_se2;
    im_msac_avg(EE,:,:) = im_msac_kern2;
    im_msac_sem(EE,:,:) = im_msac_se2;

    sim_sac_avg(EE,:,:) = sim_stim_kern;
    sim_sac_sem(EE,:,:) = sim_stim_se;
    sim_osac_avg(EE,:,:) = sim_onstim_kern2;  
    sim_osac_sem(EE,:,:) = sim_onstim_se2;
    sim_rsac_avg(EE,:,:) = sim_offstim_kern2;
    sim_rsac_sem(EE,:,:) = sim_offstim_se2;
    sim_msac_avg(EE,:,:) = sim_msac_kern2;
    sim_msac_sem(EE,:,:) = sim_msac_se2;
end

avg_gray_osac = squeeze(mean(gray_osac_avg(1:3,:,:),1));
avg_gray_rsac = squeeze(mean(gray_rsac_avg(1:3,:,:),1));
avg_gray_msac = squeeze(mean(gray_msac_avg(1:3,:,:),1));

avg_im_osac = squeeze(mean(im_osac_avg(1:3,:,:),1));
avg_im_rsac = squeeze(mean(im_rsac_avg(1:3,:,:),1));
avg_im_msac = squeeze(mean(im_msac_avg(1:3,:,:),1));

avg_sim_osac = squeeze(mean(sim_osac_avg(1:3,:,:),1));
avg_sim_rsac = squeeze(mean(sim_rsac_avg(1:3,:,:),1));
avg_sim_msac = squeeze(mean(sim_msac_avg(1:3,:,:),1));

avg_gray_sac = 0.5*avg_gray_osac + 0.5*avg_gray_rsac;
avg_im_sac = 0.5*avg_im_osac + 0.5*avg_im_rsac;

%%
load ./event_trg_avg_lfps_all_gray2.mat
gray_out_avgs = outsac_trig_avgs;
gray_ret_avgs = retsac_trig_avgs;
gray_msac_avgs = msac_trig_avgs;
gray_on_avgs = on_trig_avgs;
gray_off_avgs = off_trig_avgs;
gray_out_sems = outsac_trig_sem;
gray_ret_sems = retsac_trig_sem;
gray_msac_sems = msac_trig_sem;
gray_on_sems = on_trig_sem;
gray_off_sems = off_trig_sem;
gray_nlfps = size(gray_out_avgs,3);

load ./event_trg_avg_lfps_all_im2.mat
im_out_avgs = outsac_trig_avgs;
im_ret_avgs = retsac_trig_avgs;
im_msac_avgs = msac_trig_avgs;
im_on_avgs = on_trig_avgs;
im_off_avgs = off_trig_avgs;
im_out_sems = outsac_trig_sem;
im_ret_sems = retsac_trig_sem;
im_msac_sems = msac_trig_sem;
im_on_sems = on_trig_sem;
im_off_sems = off_trig_sem;
im_nlfps = size(im_out_avgs,3);

% load ./event_trg_avg_lfps_all_sim2.mat
% sim_out_avgs = on_trig_avgs;
% sim_ret_avgs = off_trig_avgs;
% sim_msac_avgs = msac_trig_avgs;
% sim_out_sems = on_trig_sem;
% sim_ret_sems = off_trig_sem;
% sim_msac_sems = msac_trig_sem;
% sim_nlfps = size(sim_out_avgs,3);

gray_sac_avgs = 0.5*gray_out_avgs + 0.5*gray_ret_avgs;
im_sac_avgs = 0.5*im_out_avgs + 0.5*im_ret_avgs;
gray_sac_sems = sqrt(gray_out_sems.^2 + gray_ret_sems.^2);
im_sac_sems = sqrt(im_out_sems.^2 + im_ret_sems.^2);

%%
load ./unit_ori_models
SDIM = 9;
close all
for cc = 1:96
    cc
    f1 = figure(1);
   subplot(2,4,1);hold on
   shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(1,cc,:)),squeeze(gray_sac_sem(1,cc,:)),{'b.-'},1);
   shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(1,cc,:)),squeeze(im_sac_sem(1,cc,:)),{'r.-'},1);
   axis tight
   xl = xlim();
   xlim([-0.2 0.4]);
   line(xl,[0 0],'color','k')
   subplot(2,4,2);hold on
   shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(2,cc,:)),squeeze(gray_sac_sem(2,cc,:)),{'b.-'},1);
   shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(2,cc,:)),squeeze(im_sac_sem(2,cc,:)),{'r.-'},1);
   axis tight
   xl = xlim();
   xlim([-0.2 0.4]);
   line(xl,[0 0],'color','k')
   subplot(2,4,3);hold on
   shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(3,cc,:)),squeeze(gray_sac_sem(3,cc,:)),{'b.-'},1);
   shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(3,cc,:)),squeeze(im_sac_sem(3,cc,:)),{'r.-'},1);
   axis tight
   xl = xlim();
   xlim([-0.2 0.4]);
   line(xl,[0 0],'color','k')
   subplot(2,4,4);hold on
   shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(4,cc,:)),squeeze(gray_sac_sem(4,cc,:)),{'b.-'},1);
   shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(4,cc,:)),squeeze(im_sac_sem(4,cc,:)),{'r.-'},1);
   axis tight
   xl = xlim();
   xlim([-0.2 0.4]);
   line(xl,[0 0],'color','k')
   
   subplot(2,4,5);hold on
   shadedErrorBar(lags/Fsd,squeeze(gray_sac_avgs(1,:,cc)),squeeze(gray_sac_sems(1,:,cc)),{'b-'},1);
   shadedErrorBar(lags/Fsd,squeeze(im_sac_avgs(1,:,cc)),squeeze(im_sac_sems(1,:,cc)),{'r-'},1);
   plot(lags/Fsd,gray_out_avgs(1,:,cc),'b')
   plot(lags/Fsd,gray_ret_avgs(1,:,cc),'b--')
   plot(lags/Fsd,im_out_avgs(1,:,cc),'r')
   plot(lags/Fsd,im_ret_avgs(1,:,cc),'r--')   
   axis tight
   xlim([-0.2 0.4]);
   xl = xlim();
   line(xl,[0 0],'color','k')
   
   subplot(2,4,6);hold on
   shadedErrorBar(lags/Fsd,squeeze(gray_sac_avgs(2,:,cc)),squeeze(gray_sac_sems(2,:,cc)),{'b-'},1);
   shadedErrorBar(lags/Fsd,squeeze(im_sac_avgs(2,:,cc)),squeeze(im_sac_sems(2,:,cc)),{'r-'},1);
   plot(lags/Fsd,gray_out_avgs(2,:,cc),'b')
   plot(lags/Fsd,gray_ret_avgs(2,:,cc),'b--')
   plot(lags/Fsd,im_out_avgs(2,:,cc),'r')
   plot(lags/Fsd,im_ret_avgs(2,:,cc),'r--')
   axis tight
   xlim([-0.2 0.4]);
   xl = xlim();
   line(xl,[0 0],'color','k')
   
   subplot(2,4,7);hold on
   shadedErrorBar(lags/Fsd,squeeze(gray_sac_avgs(3,:,cc)),squeeze(gray_sac_sems(3,:,cc)),{'b-'},1);
   shadedErrorBar(lags/Fsd,squeeze(im_sac_avgs(3,:,cc)),squeeze(im_sac_sems(3,:,cc)),{'r-'},1);
   plot(lags/Fsd,gray_out_avgs(3,:,cc),'b')
   plot(lags/Fsd,gray_ret_avgs(3,:,cc),'b--')
   plot(lags/Fsd,im_out_avgs(3,:,cc),'r')
   plot(lags/Fsd,im_ret_avgs(3,:,cc),'r--')
   axis tight
   xlim([-0.2 0.4]);
  xl = xlim();
   line(xl,[0 0],'color','k')
   
   subplot(2,4,8);hold on
   shadedErrorBar(lags/Fsd,squeeze(gray_sac_avgs(4,:,cc)),squeeze(gray_sac_sems(4,:,cc)),{'b-'},1);
   shadedErrorBar(lags/Fsd,squeeze(im_sac_avgs(4,:,cc)),squeeze(im_sac_sems(4,:,cc)),{'r-'},1);
   plot(lags/Fsd,gray_out_avgs(4,:,cc),'b')
   plot(lags/Fsd,gray_ret_avgs(4,:,cc),'b--')
   plot(lags/Fsd,im_out_avgs(4,:,cc),'r')
   plot(lags/Fsd,im_ret_avgs(4,:,cc),'r--')
   axis tight
   xlim([-0.2 0.4]);
  xl = xlim();
   line(xl,[0 0],'color','k')

   
   small = 0;
   big = 0;
   for i = 1:4
       subplot(2,4,i)
       cur_yl = ylim();
       if cur_yl(1) < small
           small = cur_yl(1);
       end
       if cur_yl(2) > big
           big = cur_yl(2);
       end
   end
   for i = 1:4
       subplot(2,4,i)
       ylim([small big])
   end
   
   set(f1,'Position',[100 650 1350 550])
   
   f2 = figure(2);
   cur_k = get_k_mat(ori_fit(cc));
   cur_k = reshape(cur_k,flen,SDIM);
   imagesc(un_bar_oris,1:flen,cur_k);
   mm = max(abs(cur_k(:)));
   caxis([-0.9*mm 0.9*mm])
   set(f2,'Position',[1350 650 500 400])

   
   f3 = figure(3);
     
   subplot(2,4,1);hold on
   shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(1,cc,:)),squeeze(gray_msac_sem(1,cc,:)),{'b.-'},1);
   shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(1,cc,:)),squeeze(im_msac_sem(1,cc,:)),{'r.-'},1);
   axis tight
   xl = xlim();
   line(xl,[0 0],'color','k')
   xlim([-0.2 0.4]);
   subplot(2,4,2);hold on
   shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(2,cc,:)),squeeze(gray_msac_sem(2,cc,:)),{'b.-'},1);
   shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(2,cc,:)),squeeze(im_msac_sem(2,cc,:)),{'r.-'},1);
   axis tight
   xl = xlim();
   line(xl,[0 0],'color','k')
   xlim([-0.2 0.4]);
   subplot(2,4,3);hold on
   shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(3,cc,:)),squeeze(gray_msac_sem(3,cc,:)),{'b.-'},1);
   shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(3,cc,:)),squeeze(im_msac_sem(3,cc,:)),{'r.-'},1);
   axis tight
   xl = xlim();
   line(xl,[0 0],'color','k')
   xlim([-0.2 0.4]);
   subplot(2,4,4);hold on
   shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(4,cc,:)),squeeze(gray_msac_sem(4,cc,:)),{'b.-'},1);
   shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(4,cc,:)),squeeze(im_msac_sem(4,cc,:)),{'r.-'},1);
   axis tight
   xl = xlim();
   line(xl,[0 0],'color','k')
   xlim([-0.2 0.4]);

   
   subplot(2,4,5);hold on
   shadedErrorBar(lags/Fsd,squeeze(gray_msac_avgs(1,:,cc)),squeeze(gray_msac_sems(1,:,cc)),{'b-'},1);
   shadedErrorBar(lags/Fsd,squeeze(im_msac_avgs(1,:,cc)),squeeze(im_msac_sems(1,:,cc)),{'r-'},1);
   axis tight
   xl = xlim();
   xlim([-0.2 0.4]);
   line(xl,[0 0],'color','k')
   subplot(2,4,6);hold on
   shadedErrorBar(lags/Fsd,squeeze(gray_msac_avgs(2,:,cc)),squeeze(gray_msac_sems(2,:,cc)),{'b-'},1);
   shadedErrorBar(lags/Fsd,squeeze(im_msac_avgs(2,:,cc)),squeeze(im_msac_sems(2,:,cc)),{'r-'},1);
   axis tight
   xl = xlim();
   xlim([-0.2 0.4]);
   line(xl,[0 0],'color','k')
   subplot(2,4,7);hold on
   shadedErrorBar(lags/Fsd,squeeze(gray_msac_avgs(3,:,cc)),squeeze(gray_msac_sems(3,:,cc)),{'b-'},1);
   shadedErrorBar(lags/Fsd,squeeze(im_msac_avgs(3,:,cc)),squeeze(im_msac_sems(3,:,cc)),{'r-'},1);
   axis tight
   xl = xlim();
   xlim([-0.2 0.4]);
   line(xl,[0 0],'color','k')
   subplot(2,4,8);hold on
   shadedErrorBar(lags/Fsd,squeeze(gray_msac_avgs(4,:,cc)),squeeze(gray_msac_sems(4,:,cc)),{'b-'},1);
   shadedErrorBar(lags/Fsd,squeeze(im_msac_avgs(4,:,cc)),squeeze(im_msac_sems(4,:,cc)),{'r-'},1);
   axis tight
   xl = xlim();
   xlim([-0.2 0.4]);
   line(xl,[0 0],'color','k')

      small = 0;
   big = 0;
   for i = 1:4
       subplot(2,4,i)
       cur_yl = ylim();
       if cur_yl(1) < small
           small = cur_yl(1);
       end
       if cur_yl(2) > big
           big = cur_yl(2);
       end
   end
   for i = 1:4
       subplot(2,4,i)
       ylim([small big])
   end

      set(f3,'Position',[100 0 1350 550])

   pause
   close all
end

%%
cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

cd ~/Data/bruce/G081

%fit smoothed retinotopic surface
interp_x = nan(96,1);
interp_y = nan(96,1);
tempinds = zeros(10,10);
for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
    tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
    tempinds(Y_pos(i),X_pos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');
used_inds = find(weights == 1);
tempinds = tempinds(used_inds);
interp_x(tempinds) = xpos_interp(used_inds);
interp_y(tempinds) = ypos_interp(used_inds);
xi = linspace(min(interp_x),max(interp_x),50);
yi = linspace(min(interp_y),max(interp_y),50);
[Xi,Yi] = meshgrid(xi,yi);

id_mat = nan(10,10);
for i = 1:10
    for j = 1:10
        cur = find(X_pos==j&Y_pos==i);
        if ~isempty(cur)
            id_mat(i,j) = cur;
        end
    end
end
use_ids = find(~isnan(id_mat));

%%
% [~,ord] = sort(Y_pos);
% [~,ord] = sort(X_pos);
% close all
% 
el_fdist = sqrt(interp_x.^2 + interp_y.^2);
[temp,ord] = sort(el_fdist);

% cent_point = [mean(interp_x) mean(interp_y)];
% el_edgedist = sqrt((interp_x-cent_point(1)).^2 + (interp_y-mean(cent_point(2))).^2);
% [temp,ord] = sort(el_edgedist);


figure
for i = 1:4
    subplot(2,4,i)
    cur_set = squeeze(im_msac_avgs(i,:,:))';
%     cur_set = squeeze(gray_out_avgs(i,:,:))';
%     cur_set = squeeze(sim_out_avgs(i,:,:))';
%     cur_set = squeeze(gray_on_avgs(i,:,:))';
    imagesc(lags/Fsd,1:96,cur_set(ord,:)); set(gca,'ydir','normal');
    caxis([-1.5 1.5]*1e-5)
    subplot(2,4,i+4)
    cur_set = squeeze(gray_msac_avgs(i,:,:))';
%     cur_set = squeeze(gray_ret_avgs(i,:,:))';
%     cur_set = squeeze(gray_ret_avgs(i,:,:))';
%     cur_set = squeeze(sim_ret_avgs(i,:,:))';
%     cur_set = squeeze(gray_off_avgs(i,:,:))';
    imagesc(lags/Fsd,1:96,cur_set(ord,:)); set(gca,'ydir','normal');
    caxis([-1.5 1.5]*1e-5)
end

% close all
% figure
% subplot(2,1,1)
% % cur_set = avg_im_msac;
% cur_set = avg_im_rsac;
% imagesc(tent_centers*dt,1:96,cur_set(ord,:)); set(gca,'ydir','normal');
% subplot(2,1,2)
% % cur_set = avg_gray_msac;
% cur_set = avg_gray_rsac;
% imagesc(tent_centers*dt,1:96,cur_set(ord,:)); set(gca,'ydir','normal');
% 
%%
ll = find(lags/Fsd > 0.11,1,'first');
for i = 1:4
    subplot(2,4,i)
    cur_set = squeeze(sim_out_avgs(i,ll,:));
%     cur_mat = nan(10,10);
%     cur_mat(use_ids) = cur_set(id_mat(use_ids));
%     imagesc(cur_mat); set(gca,'ydir','normal');

    F = TriScatteredInterp(interp_x,interp_y,cur_set);
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
    caxis([-1.5 1.5]*1e-5)
    
    
    subplot(2,4,i+4)
    cur_set = squeeze(gray_msac_avgs(i,ll,:));
%     cur_mat = nan(10,10);
%     cur_mat(use_ids) = cur_set(id_mat(use_ids));
%     imagesc(cur_mat); set(gca,'ydir','normal');
    F = TriScatteredInterp(interp_x,interp_y,cur_set);
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
    caxis([-1.5 1.5]*1e-5)
end

%%
close all
for i = 1:260
cur_set = squeeze(gray_msac_avgs(1,i,:));
plot(el_fdist,cur_set,'o');hold on
cur_set = squeeze(im_msac_avgs(1,i,:));
plot(el_fdist,cur_set,'ro')
title(sprintf('Time %.3f', lags(i)/Fsd));
ylim([-2 2]*1e-5);xlim([0.35 0.8]);
pause(0.1)
clf
end

%%
avg_im_msac_lfps = squeeze(mean(im_msac_avgs));
avg_gray_msac_lfps = squeeze(mean(gray_msac_avgs));
avg_sim_sac_lfps =  0.5*squeeze(mean(sim_out_avgs)) + 0.5*squeeze(mean(sim_ret_avgs));
avg_gray_sac_lfps =  0.5*squeeze(mean(gray_out_avgs)) + 0.5*squeeze(mean(gray_ret_avgs));
avg_im_sac_lfps =  0.5*squeeze(mean(im_out_avgs)) + 0.5*squeeze(mean(im_ret_avgs));
% avg_sim_sac_lfps =  0.5*squeeze(mean(im_on_avgs)) + 0.5*squeeze(mean(im_off_avgs));
% avg_sim_sac_lfps =  0.5*squeeze(mean(im_out_avgs)) + 0.5*squeeze(mean(im_ret_avgs));
% avg_sim_sac_lfps =  squeeze((im_off_avgs(4,:,:)));
close all
for i = 1:2:260
    subplot(2,1,1)
%     cur_set = squeeze(sim_out_avgs(2,i,:));
    cur_set = squeeze(avg_gray_sac_lfps(i,:))';
%     cur_set = squeeze(avg_gray_msac_lfps(i,:))';
    F = TriScatteredInterp(interp_x,interp_y,cur_set);
    Vq = F(Xi,Yi);
    pcolor(xi,yi,Vq); shading flat; colorbar;
    caxis([-1.5 1.5]*1e-5)
    title(sprintf('Time %.3f', lags(i)/Fsd));
    
    subplot(2,1,2)
%     cur_set = squeeze(gray_out_avgs(2,i,:));
    cur_set = squeeze(avg_im_sac_lfps(i,:))';
    F = TriScatteredInterp(interp_x,interp_y,cur_set);
    Vq = F(Xi,Yi);
    pcolor(xi,yi,Vq); shading flat; colorbar;
    caxis([-1.5 1.5]*1e-5)
%     title(sprintf('Time %.3f', lags(i)/Fsd));

%     subplot(3,1,3)
% %     cur_set = squeeze(gray_out_avgs(2,i,:));
%     cur_set = squeeze(avg_im_sac_lfps(i,:))';
%     F = TriScatteredInterp(interp_x,interp_y,cur_set);
%     Vq = F(Xi,Yi);
%     pcolor(xi,yi,Vq); shading flat; colorbar;
%     caxis([-1.5 1.5]*1e-5)

pause(0.1)
    clf
end
