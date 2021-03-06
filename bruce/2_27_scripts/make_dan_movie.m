
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
dsfrac = 0.5;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = Nyp/dsfrac;
Nx = Nxp/dsfrac;
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [4 5; -3.5 -2.5]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [-0.5 0.5; -0.5 0.5]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

% desired temporal resolution
stimres = 0.025; %in s
deltaT = stimres*1; % frame rate for movie

% some options
SHOWMOVIE = 1; % show movie if 1
SAVEDATA = 0;  % save data to disk if 1

% block number
% blockid = 4;
blockid = 1
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
% for i=1:Nimage
%%
close all
i = 11;%10

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

% IMAGE = IMAGEorg(:,129:1152);
IMAGE = IMAGEorg;
IMAGE = imresize(IMAGE,2); %image downsampling
[Nyd Nxd] = size(IMAGE); %down-sampled image size

%if you want to reflect the image about x or y axes
IMAGE = flipud(IMAGE);
% imagesc(IMAGE);colormap(gray);set(gca,'ydir','normal');
%%

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
lEyelow = interp1(eyets(eyeindx), lEye(eyeindx,:), Tticks);
%     lEyelow = interp1(eyets(eyeindx), rEye(eyeindx,:), Tticks);

%allow for extrapolation of eye signal to first or last sample if
%interpolation failed
if any(isnan(lEyelow(end,:)))
    lEyelow(end,:) = lEyelow(end-1,:);
end
if any(isnan(lEyelow(1,:)))
    lEyelow(1,:) = lEyelow(2,:);
end

% convert eye signal to pix
eyepix = lEyelow./Pix2Deg/dsfrac;

%%
% shift the image based on eye signal
NT = size(lEyelow,1);
%     if SAVEDATA
STIMrec = zeros([NT size(IMAGE)]);
%     end
cd ~/dan_movie2/
gaze_buff = [];
buff_size = NT;
for frame=1:NT
    gaze_buff = [lEyelow(frame,:); gaze_buff];
    if size(gaze_buff,1) > buff_size
        gaze_buff(end,:) = [];
    end
    % shift image based on current eye position
    STstim = dist_shift2d(IMAGE,-eyepix(frame,1),2,0);
    STstim = dist_shift2d(STstim,-eyepix(frame,2),1,0);
    STstim_patch = STstim(ypatch_inds,xpatch_inds);
    
    if SHOWMOVIE
        %             [Ny,Nx] = size(IMAGE);
        figure
        colormap gray;
        set(gca,'fontsize',14)
        %static image
        subplot(2,1,1);
        imagesc(xax,yax,IMAGE); title('Static Image','fontsize',16); set(gca,'ydir','normal')
        hold on
%         axis square
        xlabel('Horizontal Position (degrees)','fontsize',14)
        ylabel('Vertical Position (degrees)','fontsize',14)
%         scatter(gaze_buff(:,1),gaze_buff(:,2),[buff_size:-1:(buff_size-size(gaze_buff,1)+1)]'*2,'ro')
        plot(gaze_buff(:,1),gaze_buff(:,2),'r-','linewidth',1)
        le = RF_patch(1,1)+lEyelow(frame,1);
        re = le + RF_patch_width(1);
        be = RF_patch(2,1)+lEyelow(frame,2);
        ue = be + RF_patch_width(2);
        line([le re],(ue+be)/2+[0 0],'color','w','linewidth',1)
        line((le+re)/2+[0 0],[be ue],'color','w','linewidth',1)
         r= rectangle('Position',[RF_patch(1,1)+lEyelow(frame,1) RF_patch(2,1)+lEyelow(frame,2) RF_patch_width(1) RF_patch_width(2)],...
            'EdgeColor','r','LineWidth',4);
       
        %RF patch stimulus
        subplot(2,1,2);
        %             imagesc(xax,yax,STstim); title('Real Stimulus'); set(gca,'ydir','normal')
        %             ylim(RF_patch(2,:)); xlim(RF_patch(1,:));
        imagesc(xax(xpatch_inds),yax(ypatch_inds),STstim_patch); title('Image on retina','fontsize',16); set(gca,'ydir','normal')
        xl = xlim();
        yl = ylim();
        line(xl,[0 0],'color','w','linewidth',1)
        line([0 0],yl,'color','w','linewidth',1)
        axis square
        set(gca,'xtick',[-0.5:0.1:0.5])
        ylim(RF_patch(2,:)); xlim(RF_patch(1,:));
        xlabel('Horizontal Position (degrees)','fontsize',14)
        ylabel('Vertical Position (degrees)','fontsize',14)

        fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [7 12]);
        fname = sprintf('dan_movie_%d',333+frame);
        print('-dtiff',fname);close 
%         allpause(.1); clf
    end
end

% end
