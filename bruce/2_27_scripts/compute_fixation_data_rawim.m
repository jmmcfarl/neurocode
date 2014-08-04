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
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

% desired temporal resolution
stimres = 0.025; %in s
deltaT = stimres*2; % frame rate for movie

gaussfilt = fspecial('gaussian',5,0.25);

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v4
cd ~/Data/bruce/2_27_12/ExptA
%% IF YOU WANT DISPARITY BASED GAIN CORRECTION
% actual_disp = all_sac_locs_right - all_sac_locs_left;
% med_disp = median(actual_disp);
% all_sac_locs_right = bsxfun(@minus,all_sac_locs_right,med_disp);
% actual_disp = all_sac_locs_right - all_sac_locs_left;
%
% K0 = [1 0 0 1 1 0 0 1];
% options.Display = 'none';
% scale_pen = 1;
% [Kbest] = fminunc(@(K) disparity_error(K,all_sac_locs_right,all_sac_locs_left,scale_pen), K0,options);
% clear cur_r cur_l
% cur_r(:,1) = Kbest(1)*all_sac_locs_right(:,1) + Kbest(2)*all_sac_locs_right(:,2);
% cur_r(:,2) = Kbest(3)*all_sac_locs_right(:,1) + Kbest(4)*all_sac_locs_right(:,2);
%
% cur_l(:,1) = Kbest(5)*all_sac_locs_left(:,1) + Kbest(6)*all_sac_locs_left(:,2);
% cur_l(:,2) = Kbest(7)*all_sac_locs_left(:,1) + Kbest(8)*all_sac_locs_left(:,2);
%
% cur_disp = cur_r - cur_l;
%
% [rbefore,b] = corr(actual_disp, all_sac_locs_right);
% [lbefore,b] = corr(actual_disp, all_sac_locs_left);
%
% [rafter,b] = corr(cur_disp, cur_r);
% [lafter,b] = corr(cur_disp, cur_l);
%
% fprintf('Before: L=%.3f  R=%.3f\n',sum(abs(lbefore(:))),sum(abs(rbefore(:))));
% fprintf('After: L=%.3f  R=%.3f\n',sum(abs(lafter(:))),sum(abs(rafter(:))));
%
% all_sac_locs_right = cur_r;
% all_sac_locs_left = cur_l;

%%
% [~,hsort] = sort(all_sac_locs_left(:,1));
% figure
% plot(all_sac_locs_left(hsort,1),actual_disp(hsort,1),'.')
% hold on
% plot(all_sac_locs_left(hsort,1),smooth(actual_disp(hsort,1),20,'rlowess'),'r','linewidth',2)
% %
% [~,hsort] = sort(all_sac_locs_right(:,1));
% figure
% plot(all_sac_locs_right(hsort,1),actual_disp(hsort,1),'.')
% hold on
% plot(all_sac_locs_right(hsort,1),smooth(actual_disp(hsort,1),20,'rlowess'),'r','linewidth',2)
%
% [~,hsort] = sort(cur_l(:,1));
% figure
% plot(cur_l(hsort,1),cur_disp(hsort,1),'.')
% hold on
% plot(cur_l(hsort,1),smooth(cur_disp(hsort,1),20,'rlowess'),'k','linewidth',2)
% %
% [~,hsort] = sort(cur_r(:,1));
% figure
% plot(cur_r(hsort,1),cur_disp(hsort,1),'.')
% hold on
% plot(cur_r(hsort,1),smooth(cur_disp(hsort,1),20,'rlowess'),'k','linewidth',2)

%%
used_fixs = find(max(abs(all_sac_locs_right_nt),[],2) < 5 & max(abs(all_sac_locs_left_nt),[],2) < 5);
n_fixs = length(used_fixs);

X_left = zeros(n_fixs,length(ypatch_inds),length(xpatch_inds));
X_right = zeros(n_fixs,length(ypatch_inds),length(xpatch_inds));
X_avg = zeros(n_fixs,length(ypatch_inds),length(xpatch_inds));

% all_sac_locs_left_pix = round(all_sac_locs_left/Pix2Deg/dsfrac);
% all_sac_locs_right_pix = round(all_sac_locs_right/Pix2Deg/dsfrac);
all_sac_locs_left_pix = round(all_sac_locs_left_nt/Pix2Deg/dsfrac);
all_sac_locs_right_pix = round(all_sac_locs_right_nt/Pix2Deg/dsfrac);
all_sac_locs_avg_pix = 0.5*all_sac_locs_left_pix + 0.5*all_sac_locs_right_pix;
for blockid = 1:5
    
    stimtime = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    Nimage = length(stimID);
    
    % load images
    for i=1:Nimage
        
        cur_fixs = find(all_stim_nums(used_fixs) == i & blockids(used_fixs) == blockid);
        if ~isempty(cur_fixs)
            fprintf('Image %d of %d\n',i,Nimage);
            clear filename 
            if stimID(i)<10
                filename = ['ExptA0000' num2str(stimID(i)) '.png'];
            elseif stimID(i)>=10 && stimID(i)<100
                filename = ['ExptA000' num2str(stimID(i)) '.png'];
            elseif stimID(i) >= 100
                filename = ['ExptA00' num2str(stimID(i)) '.png'];
            end
            
            IMAGEorg = imread(filename);
            IMAGEorg = double(IMAGEorg); % convert to double format
            IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
            
            [Ny Nx] = size(IMAGEorg);
            
            IMAGE = IMAGEorg;
            %             IMAGE = imfilter(IMAGEorg,gaussfilt,'replicate');%slight smoothing filter (anti-aliasing)
            IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
            [Nyd Nxd] = size(IMAGE); %down-sampled image size
            
            %if you want to reflect the image about x or y axes
            IMAGE = flipud(IMAGE);
            %     IMAGE = fliplr(IMAGE);
            
            %zero pad image
            nzpad = 200;
            IMAGEp = cat(1,nan(nzpad,size(IMAGE,2)),IMAGE);
            IMAGEp = cat(1,IMAGEp,nan(nzpad,size(IMAGEp,2)));
            IMAGEp = cat(2,nan(size(IMAGEp,1),nzpad),IMAGEp);
            IMAGEp = cat(2,IMAGEp,nan(size(IMAGEp,1),nzpad));
            
            for j = 1:length(cur_fixs)
                ypatch_inds_adj = round(ypatch_inds + all_sac_locs_left_pix(used_fixs(cur_fixs(j)),2) + nzpad);
                xpatch_inds_adj = round(xpatch_inds + all_sac_locs_left_pix(used_fixs(cur_fixs(j)),1) + nzpad);
                STstim_patch = IMAGEp(ypatch_inds_adj,xpatch_inds_adj);
                
                X_left(cur_fixs(j),:,:) = STstim_patch - mean(STstim_patch(:));
                
                ypatch_inds_adj = round(ypatch_inds + all_sac_locs_right_pix(used_fixs(cur_fixs(j)),2) + nzpad);
                xpatch_inds_adj = round(xpatch_inds + all_sac_locs_right_pix(used_fixs(cur_fixs(j)),1) + nzpad);
                STstim_patch = IMAGEp(ypatch_inds_adj,xpatch_inds_adj);
                
                X_right(cur_fixs(j),:,:) = STstim_patch - mean(STstim_patch(:));
                
                ypatch_inds_adj = round(ypatch_inds + all_sac_locs_avg_pix(used_fixs(cur_fixs(j)),2) + nzpad);
                xpatch_inds_adj = round(xpatch_inds + all_sac_locs_avg_pix(used_fixs(cur_fixs(j)),1) + nzpad);
                STstim_patch = IMAGEp(ypatch_inds_adj,xpatch_inds_adj);
                
                X_avg(cur_fixs(j),:,:) = STstim_patch - mean(STstim_patch(:));
            end
        end
    end
end

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save fixation_image_patches_d4_nt_origcal_2 X_left X_right X_avg used_fixs