clear all
% close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;


Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

% desired temporal resolution
stimres = 0.025; %in s
deltaT = stimres*2; % frame rate for movie

% gaussfilt = fspecial('gaussian',5,0.25);

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
%%
complexity_index_fin = log10(abs(gabor_params_fin(:,7))./sqrt(sum(gabor_params_fin(:,8:9).^2,2)));
complexity_index_ini = log10(abs(gabor_params{1}(:,7))./sqrt(sum(gabor_params{1}(:,8:9).^2,2)));
muaids = [1 2 4 5 6 11 12 13 14];
cellids = [1 2 3 4 5 6 7 8 10];
cd /Users/James/James_scripts/bruce/gabor_models
    close all
for n = 1:18;
    
    cur_mask_init1 = get_pgabor_mask(gabor_params{1}(n,1:6),0,[SDIM SDIM]);
    cur_mask_init2 = get_pgabor_mask(gabor_params{1}(n,1:6),pi/2,[SDIM SDIM]);
    cur_mask_init = gabor_params{1}(n,8)*cur_mask_init1 + gabor_params{1}(n,9)*cur_mask_init2;
    
    cur_mask_fin1 = get_pgabor_mask(gabor_params_fin(n,1:6),0,[SDIM SDIM]);
    cur_mask_fin2 = get_pgabor_mask(gabor_params_fin(n,1:6),pi/2,[SDIM SDIM]);
    cur_mask_rev = gabor_params{end}(n,8)*cur_mask_fin1 + gabor_params{end}(n,9)*cur_mask_fin2;
    f1 = figure;
    subplot(2,2,1)
    imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_mask_init);
    set(gca,'ydir','normal');
    xlabel('Xposition (degrees)','fontsize',14)
    ylabel('Yposition (degrees)','fontsize',14)
    title('Initial fit','fontsize',16)
    
    subplot(2,2,2)
    imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_mask_rev)
    set(gca,'ydir','normal');
     xlabel('Xposition (degrees)','fontsize',14)
    ylabel('Yposition (degrees)','fontsize',14)
    title('Final fit','fontsize',16)
   
        subplot(2,2,[3 4])
    if n <= length(cellids)
        plot(spk_cnts(:,cellids(n)))
        axis tight
        title(sprintf('CI: %.2f',complexity_index_fin(n)));
        
        fname = sprintf('Cell%d_Modelfits',cellids(n));
        
    else
        cur_n = n - length(cellids);
        plot(mua_cnts(:,muaids(cur_n)))
        axis tight
        title(sprintf('CI: %.2f',complexity_index_fin(n)));
        fname = sprintf('MUA%d_Modelfits',muaids(cur_n));
    end
        xlabel('Fixation number','fontsize',14)
        ylabel('Spk cnt','fontsize',14)
    
    fp = fillPage(f1, 'margins', [0 0 0 0], 'papersize', [10 10]);
    print('-dpdf','-painters',fname);close all
end

% pref_ori = gabor_params_fin(n,3)*180/pi;
% pref_freq = 1/gabor_params_fin(n,4);
% phases = linspace(0,360,50);
% for i = 1:length(phases)
%     [ipat] = hartley3p(pref_ori,pref_freq*RF_patch_width(1)*Fsd,phases(i),SDIM);
%     m1_out = dot(cur_mask_fin1(:),ipat(:));
%     m2_out = dot(cur_mask_fin2(:),ipat(:));
%     e_out = sqrt(m1_out.^2+m2_out.^2);
%     gout = gabor_params_fin(n,7)*e_out + gabor_params_fin(n,8)*m1_out + gabor_params_fin(n,9)*m2_out;
%     rpred(i) = log(1+exp(gout+gabor_params_fin(n,10)));
% end
% roffset = log(1+exp(gabor_params_fin(n,10)));
% subplot(2,2,3)
% plot(phases,rpred-roffset)
% axis tight
