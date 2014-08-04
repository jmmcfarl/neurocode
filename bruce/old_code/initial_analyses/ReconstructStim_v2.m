% reconstruct the stimulus
% created by Yuwei Cui, Feb 28, 2012
% edited by James McFarland

clear all close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/

%%
load Blocks.mat

% original image resolution
Pix2Deg = 1.1279 / 60;

% down-sampling fraction for image
dsfrac = 4;

% desired temporal resolution
stimres = 0.05; %in s

% some options
SHOWMOVIE = 1; % show movie if 1
SAVEDATA = 0;  % save data to disk if 1

% block number
blockid = 1;

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

gaussfilt = fspecial('gaussian',5,0.5);

ov_stimRec = [];
ov_eyepos = [];

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
%     IMAGE = flipud(IMAGE);
%     IMAGE = fliplr(IMAGE);
    
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
%     lEyelow = interp1(eyets(eyeindx), lEye(eyeindx,:), Tticks);
    lEyelow = interp1(eyets(eyeindx), rEye(eyeindx,:), Tticks);
    
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
    deltaT = stimres*2; % frame rate for movie
%     if SAVEDATA
        STIMrec = zeros([NT size(IMAGE)]);
%     end
        
    for frame=1:NT
        % shift image based on current eye position
        STstim = dist_shift2d(IMAGE,-eyepix(frame,1),2,0);
        STstim = dist_shift2d(STstim,eyepix(frame,2),1,0);
        if SAVEDATA
            STIMrec(frame,:,:) = STstim;
        end        
        
        if SHOWMOVIE
            [Nx,Ny] = size(IMAGE);
            figure(97);clf;
            colormap gray;
            subplot(2,2,1);
            imagesc(STstim); title('Real Stimulus');
            subplot(2,2,2);
            imagesc(IMAGE); title('Static Image');
            subplot(2,2,3);title('Eye Position');
            plot(eyepix(1:frame,1),eyepix(1:frame,2));hold on;
            plot(eyepix(frame,1),eyepix(frame,2),'rx');
            xlim([-Nx/2 Nx/2]);
            ylim([-Ny/2 Ny/2]);
            subplot(2,2,4);
            imagesc(STstim); title('Real Stimulus');
%             ylim([120 133]); xlim([150 163]);
            pause(deltaT);
            
        end
    end
    if SAVEDATA
%         save(['BLOCK' num2str(blockid) 'IMAGE' num2str(i) 'r.mat'], 'STIMrec','STIMrecH','STIMrecV');
        save(['BLOCK' num2str(blockid) 'IMAGE' num2str(i) 'r2_p5.mat'], 'STIMrec');
    end
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
% save stimrecon_t stimres stim_num recon_t n_recon_samps