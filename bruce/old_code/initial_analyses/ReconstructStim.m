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
% fracx = 4; fracy = 4;
dsfrac = 4;

% crop range (in deg) (shown as rectangular at eye position plot)
crop = [-7 7 -7 7];

% desired temporal resolution
stimres = 0.03;

% some options
SHOWMOVIE = 1; % show movie if 1
SAVEDATA = 0;  % save data to disk if 1
% block number
blockid = 1;


stimtime = Blocks{blockid}.stimtime;
stimID = Blocks{blockid}.stimids;
Nimage = length(stimID);

croppix0 = round(crop./Pix2Deg);
% croppix0(1:2) = croppix0(1:2)./fracx;
% croppix0(3:4) = croppix0(3:4)./fracy;
% load eye-movement
load( ['~/Data/bruce/2_27_12/saccades/lemM232.' num2str(50+blockid) '.em.sac.mat']);
% left eye
lEye = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
% right eye
% lEye = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
EyeEndT = Expt.Trials.End/10000; % time of first eye sample
Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
eyets = EyeStartT:Eyedt:EyeEndT;

recon_t = [];
stim_num = [];
gaussfilt = fspecial('gaussian',5,1);

[X,Y] = meshgrid(1:2:320,1:2:256);
[Xi,Yi] = meshgrid(1:320,1:256);

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
    % crop image
    
    %     croppix(3:4) = croppix0(3:4)+round(Ny/2);
    %     croppix(1:2) = croppix0(1:2)+round(Nx/2);
    %     IMAGEorg = IMAGEorg(croppix(3):croppix(4), croppix(1):croppix(2));
    
    IMAGE = imfilter(IMAGEorg,gaussfilt,'replicate');%slight smoothing filter
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    [Nyd Nxd] = size(IMAGE);
%     IMAGE = flipud(IMAGE);
%     IMAGE = fliplr(IMAGE);

    %     % down sampling image
    %     IMAGE = DownSampleImg(IMAGEorg, fracx, fracy);
    %     %     figure,imagesc(IMAGE);colormap gray
    
    onsetT = stimtime(i);
    if i < Nimage
        endT = stimtime(i+1)-stimres;
        endT = endT - 0.2; %last 200ms doesn't have any image presentation!
    else
        endT = Blocks{blockid}.blocktimes(2,end);
    end
    
    % relevant eye signal
    eyeindx = round((onsetT-EyeStartT)/Eyedt):round((endT-EyeStartT)/Eyedt);
    eyeindx(eyeindx > length(eyets) | eyeindx > size(lEye,1)) = [];
    Tticks = onsetT: stimres: endT;
    Tticks(Tticks > eyets(eyeindx(end))) = [];
    lEyelow = interp1(eyets(eyeindx), lEye(eyeindx,:), Tticks); % interpolation
    
    %allow for extrapolation of eye signal to first or last sample if
    %interpolation failed
    if any(isnan(lEyelow(end,:)))
        lEyelow(end,:) = lEyelow(end-1,:);
    end
    if any(isnan(lEyelow(1,:)))
        lEyelow(1,:) = lEyelow(2,:);
    end
    
    % to compare original eye signal and down-sampled eye signal
    %figure,plot(Tticks,lEyelow,'b');hold on;plot(eyets(eyeindx), lEye(eyeindx,:),'r');
    
    % convert eye signal to pix
    eyepix = lEyelow./Pix2Deg;
    eyepix(:,1) = eyepix(:,1)/dsfrac;
    eyepix(:,2) = eyepix(:,2)/dsfrac;
    
    % shift the image based on eye signal
    NT = size(lEyelow,1);
    deltaT = stimres*2; % frame rate for movie
    if SAVEDATA
        STIMrec = zeros([NT size(IMAGE)]);
    end
    
    %to simulate effects of uncertainty on eye movements
    %                 eyepix = eyepix + randn(size(eyepix))*0.2/(Pix2Deg*dsfrac);
    
    [cA,cH,cV,cD] = dwt2(IMAGE,'sym1','mode','sym');
    cH = imfilter(abs(cH),gaussfilt,'replicate');%slight smoothing filter
    cV = imfilter(abs(cV),gaussfilt,'replicate');%slight smoothing filter
    new_stimH = sqrt(interp2(X,Y,cH,Xi,Yi));
    new_stimV = sqrt(interp2(X,Y,cV,Xi,Yi));
    
    for frame=1:NT
        % shift image
        STstim = dist_shift2d(IMAGE,-eyepix(frame,1),2,0);
        STstim = dist_shift2d(STstim,eyepix(frame,2),1,0);
        STstimH = dist_shift2d(new_stimH,-eyepix(frame,1),2,0);
        STstimH = dist_shift2d(STstimH,eyepix(frame,2),1,0);
        STstimV = dist_shift2d(new_stimV,-eyepix(frame,1),2,0);
        STstimV = dist_shift2d(STstimV,eyepix(frame,2),1,0);
        if SAVEDATA
            STIMrec(frame,:,:) = STstim;
            STIMrecH(frame,:,:) = STstimH;
            STIMrecV(frame,:,:) = STstimV;
        end
        
        
        if SHOWMOVIE
            [Nx,Ny] = size(IMAGE);
            figure(97);clf;
            colormap gray;
            subplot(3,2,1);
            imagesc(STstim); title('Real Stimulus');
            subplot(3,2,2);
            imagesc(IMAGE); title('Static Image');
            subplot(3,2,3);title('Eye Position');
            plot(eyepix(1:frame,1),eyepix(1:frame,2));hold on;
            plot(eyepix(frame,1),eyepix(frame,2),'rx');
            xlim([-Nx/2 Nx/2]);
            ylim([-Ny/2 Ny/2]);
            subplot(3,2,4);
            imagesc(STstim); title('Real Stimulus');
%             ylim([120 133]); xlim([150 163]);
            subplot(3,2,5);
            imagesc(STstimH); title('Real Stimulus');
%             ylim([120 133]); xlim([150 163]);
            %             caxis([0 20])
caxis([1 6])            
subplot(3,2,6);
            imagesc(STstimV); title('Real Stimulus');
%             ylim([120 133]); xlim([150 163]);
            %             caxis([0 20])
caxis([1 6])            %             rectangle('Position',[croppix(1)-round(Nx/2),croppix(3)-round(Ny/2), (croppix(2)-croppix(1)), (croppix(4)-croppix(3))])
            %             axis(croppix0*1.2);
            pause(deltaT);
            
        end
    end
    if SAVEDATA
%         save(['BLOCK' num2str(blockid) 'IMAGE' num2str(i) 'r.mat'], 'STIMrec','STIMrecH','STIMrecV');
        save(['BLOCK' num2str(blockid) 'IMAGE' num2str(i) 'lfly.mat'], 'STIMrec');
    end
    
    n_recon_samps(i) = NT;
    recon_t = [recon_t stimtime(i):stimres:(stimtime(i)+n_recon_samps(i)*stimres)];
    stim_num = [stim_num i*ones(1,n_recon_samps(i))];
    
end

cd ~/Data/bruce/2_27_12/stimrecon/
save stimrecon_t stimres stim_num recon_t n_recon_samps