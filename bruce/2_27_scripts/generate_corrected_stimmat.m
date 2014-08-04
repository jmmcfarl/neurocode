
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
load fixation_data_v4

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);

%%
n_fixs = length(used_fixs);

X = zeros(n_fixs,length(ypatch_inds),length(xpatch_inds));
stim_mean = zeros(n_fixs,1);
stim_cont = zeros(n_fixs,1);

est_pos = [tot_est_X tot_est_Y];
est_pos_pix = round(est_pos/Pix2Deg/dsfrac);

for blockid = 1:5
    
    stimtime = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    Nimage = length(stimID);
    
    % load images
    for i=1:Nimage
        
        cur_fixs = find(all_stim_nums(used_fixs) == i & blockids(used_fixs) == blockid);
        if ~isempty(cur_fixs)
            clear filename
            fprintf('Image %d of %d\n',i,Nimage);
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
            nzpad = 300;
            IMAGEp = cat(1,nan(nzpad,size(IMAGE,2)),IMAGE);
            IMAGEp = cat(1,IMAGEp,nan(nzpad,size(IMAGEp,2)));
            IMAGEp = cat(2,nan(size(IMAGEp,1),nzpad),IMAGEp);
            IMAGEp = cat(2,IMAGEp,nan(size(IMAGEp,1),nzpad));
            
            for j = 1:length(cur_fixs)
                ypatch_inds_adj = round(ypatch_inds + est_pos_pix(cur_fixs(j),2) + nzpad);
                xpatch_inds_adj = round(xpatch_inds + est_pos_pix(cur_fixs(j),1) + nzpad);
                STstim_patch = IMAGEp(ypatch_inds_adj,xpatch_inds_adj);
                stim_mean(cur_fixs(j)) = mean(STstim_patch(:));
X(cur_fixs(j),:,:) = STstim_patch;
                                X(cur_fixs(j),:,:) = STstim_patch - stim_mean(cur_fixs(j));
                stim_cont(cur_fixs(j)) = std(STstim_patch(:));
                
%                 X(cur_fixs(j),:,:) = X(cur_fixs(j),:,:)/stim_cont(cur_fixs(j));
            end
        end
    end
end

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save fixation_image_patches_corrected X used_fixs est_pos* stim_cont stim_mean