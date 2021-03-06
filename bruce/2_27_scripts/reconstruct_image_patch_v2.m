
clear all close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/
addpath('~/James_scripts/bruce/');

%%
load Blocks.mat

% original image resolution
% Pix2Deg = 1.1279 / 60;
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

RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [3. 6.5; -4 -0.5]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [1.5 8; -5.5 1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

% desired temporal resolution
stimres = 0.025; %in s
deltaT = stimres*2; % frame rate for movie

% some options
SHOWMOVIE = 0; % show movie if 1
SAVEDATA = 0;  % save data to disk if 1

% block number
% blockid = 4;
for blockid = 1:4
    stimtime = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    Nimage = length(stimID);
    
    % load eye-movement
    load( ['~/Data/bruce/2_27_12/saccades/lemM232.' num2str(50+blockid) '.em.sac.mat']);
    
    % left eye
    lEye = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
    % right eye
    rEye = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
    
    recon_t = []; %time axis for reconstructed stim matrix (sec)
    stim_num = []; %corresponding vector of image indices
    
    gaussfilt = fspecial('gaussian',5,0.25);
    
    ov_stimRec = [];
    ov_eyepos = [];
    ov_im_patch = [];
    
    % load images
    for i=1:Nimage
        
        fprintf('Image %d of %d\n',i,Nimage);
        if stimID(i)<10
            filename = ['ExptA0000' num2str(stimID(i)) '.png'];
        elseif stimID(i)>=10 && stimID(i)<100
            filename = ['ExptA000' num2str(stimID(i)) '.png'];
        end
        
        IMAGEorg = imread(filename);
        IMAGEorg = double(IMAGEorg); % convert to double format
        IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
        
        [Ny Nx] = size(IMAGEorg);
        
        IMAGE = imfilter(IMAGEorg,gaussfilt,'replicate');%slight smoothing filter (anti-aliasing)
        IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
        [Nyd Nxd] = size(IMAGE); %down-sampled image size
        
        %if you want to reflect the image about x or y axes
        IMAGE = flipud(IMAGE);
            IMAGE = fliplr(IMAGE);
        
        %zero pad image
        nzpad = 200;
        IMAGEp = cat(1,nan(nzpad,size(IMAGE,2)),IMAGE);
        IMAGEp = cat(1,IMAGEp,nan(nzpad,size(IMAGEp,2)));
        IMAGEp = cat(2,nan(size(IMAGEp,1),nzpad),IMAGEp);
        IMAGEp = cat(2,IMAGEp,nan(size(IMAGEp,1),nzpad));
        
        
        onsetT = stimtime(i); %stimulus onset time (s)
        if i < Nimage %for the last image, take the end time as the block end time
            endT = stimtime(i+1)-stimres; %set end time as one time bin before the following stimulus presentation
            endT = endT - 0.2; %last 200ms doesn't have any image presentation!
        else
            endT = Blocks{blockid}.blocktimes(2,end);
        end
        
        % relevant eye signal
        eyeindx = round((onsetT-EyeStartT)/Eyedt):round((endT-EyeStartT)/Eyedt); %eye-tracking inds during current stimulus
        
        %make sure we aren't outside the range of the eye-tracking data
        eyeindx(eyeindx > length(eyets) | eyeindx > size(lEye,1)) = [];
        
        %time axis for reconstructing current image stimulus
        Tticks = onsetT: stimres: endT;
        Tticks(Tticks > eyets(eyeindx(end))) = [];
        
        %interpolate eye-tracking signal onto current time axis (linear interp)
%         lEyelow = interp1(eyets(eyeindx), lEye(eyeindx,:), Tticks);
%             lEyelow = interp1(eyets(eyeindx), rEye(eyeindx,:), Tticks);
lEyelow = interp1(eyets(eyeindx),lEye(eyeindx,:)+rEye(eyeindx,:),Tticks)/2; %average of left and right eye positions        

        %allow for extrapolation of eye signal to first or last sample if
        %interpolation failed
        if any(isnan(lEyelow(end,:)))
            lEyelow(end,:) = lEyelow(end-1,:);
        end
        if any(isnan(lEyelow(1,:)))
            lEyelow(1,:) = lEyelow(2,:);
        end
        
        ov_eyepos = [ov_eyepos; lEyelow];
        
        % convert eye signal to pix
        eyepix = lEyelow./Pix2Deg/dsfrac;
        
        % shift the image based on eye signal
        NT = size(lEyelow,1);
        %     if SAVEDATA
        STIMrec = zeros([NT size(IMAGE)]);
        %     end
        
        image_avg(i) = mean(IMAGE(:));
        patch_avg{i} = zeros(NT,1);
        for frame=1:NT
            % shift image based on current eye position
            ypatch_inds_adj = round(ypatch_inds + eyepix(frame,2) + nzpad);
            xpatch_inds_adj = round(xpatch_inds + eyepix(frame,1) + nzpad);
            STstim_patch = IMAGEp(ypatch_inds_adj,xpatch_inds_adj);
            patch_avg{i}(frame) = nanmean(STstim_patch(:));
            ov_im_patch = cat(1,ov_im_patch,reshape(STstim_patch,[1 size(STstim_patch,1) size(STstim_patch,2)]));
        end
        all_patch_avg(i) = mean(patch_avg{i});
        %     ov_stimRec = [ov_stimRec; STIMrec];
        
        n_recon_samps(i) = NT;
        %     recon_t = [recon_t stimtime(i):stimres:(stimtime(i)+(n_recon_samps(i)-1)*stimres)];
        recon_t = [recon_t Tticks];
        stim_num = [stim_num i*ones(1,n_recon_samps(i))];
        
    end
    
    recon_t_or = recon_t;
    stim_num_or = stim_num;
    n_recon_samps_or = n_recon_samps;
    cd ~/Data/bruce/2_27_12/stimrecon/
    % save used_eyesig_p5 ov_eyepos recon_t_or stim_num_or n_recon_samps_or
    % save stimrecon_t_cor stimres stim_num recon_t n_recon_samps
    sname = sprintf('image_patch_block%d_avgeye_25_dsf4_flipx',blockid);
    save(sname,'stimres','stim_num','recon_t','n_recon_samps','ov_im_patch');
    
end