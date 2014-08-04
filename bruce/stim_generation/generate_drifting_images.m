close all
clear all

save_dir = '~/Data/bruce/Expt_1_9_13_imfolder';
% save_dir = './';
blackColor = 0;
whiteColor = 255;
% cbounds = [25 75];
cbounds = [5 95];
% cbounds = [10 90];
grayColor = (blackColor+whiteColor)/2;
Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
niqf = Fs/2;
xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;
[XAX,YAX] = meshgrid(xax,yax);

% center pixel before scale transformation
XcenterPix = Nxp/2;
YcenterPix = Nyp/2;

% sac_amp = 1.5;
% sacs_per_im = 4;

N_nat_trials = 5;
ims_per_trial = 160;

noise_std = 1.5/Pix2Deg;
temp_filt_std = 8;
temp_filt_ax = -4*temp_filt_std:4*temp_filt_std;
temp_filt = exp(-temp_filt_ax.^2/(2*temp_filt_std^2));
temp_filt =temp_filt/sum(temp_filt);


%% FOR NAT IMAGE TRIALS
nat_dir = '~/James_scripts/data_processing/Images/image_set_A/'
% file_range_nat = [1:685];
% rr = range(file_range_nat)+1;
% file_order_nat = randperm(rr) + file_range_nat(1)-1;

use_files = [4 6 9 16 599];
%%
% cnt = 1;
% for cur_trial = 1:N_nat_trials
%     fprintf('Nat image trial %d of %d\n',cur_trial,N_nat_trials);
%     cur_im = 1;
%     
%     myimgfile = sprintf('%.4d.png',use_files(cnt));
%     cnt = cnt + 1;
%     
%     imdata = imread(myimgfile);
%     imdata = double(imdata); % convert to double format
%     
%     cur_im
%     imagesc(imdata);colormap(gray);%set(gca,'ydir','normal')
%     pause
% end

%%
cnt = 1;
trial_offset = 0;
cd(nat_dir)
for cur_trial = 1:N_nat_trials
    fprintf('Nat image trial %d of %d\n',cur_trial,N_nat_trials);
    cur_im = 1;
    
    myimgfile = sprintf('%.4d.png',use_files(cnt));
    cnt = cnt + 1;
    
    imdata = imread(myimgfile);
    imdata = double(imdata); % convert to double format
    
    imdataorg = imdata;
    
    %%
    x_white = randn(ims_per_trial,1)*noise_std;
    y_white = randn(ims_per_trial,1)*noise_std;
    x_filt = conv(x_white,temp_filt,'same');
    y_filt = conv(y_white,temp_filt,'same');
    x_filt = x_filt/std(x_filt)*noise_std;
    y_filt = y_filt/std(y_filt)*noise_std;
    %%
    X_inds = round([x_filt(:) y_filt(:)]);
    
    %%
    [Ny,Nx] = size(imdata);
    %% cd(save_dir)
%     colormap(gray)
    for cur_frame = 1:ims_per_trial
        shift_im = shift_mat2(imdataorg,X_inds(cur_frame,1),X_inds(cur_frame,2),1);
        imdata = shift_im;
        
%         imagesc(imdata); pause(0.025)
        
        %print raw image
        imdatapng = imdata./whiteColor;
        
        im_id = 4000000+1000*(cur_trial+trial_offset) + cur_frame;
        cur_fname = sprintf('IM%.7d',im_id);
        imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.pgm'),'pgm');
        cur_im = cur_im + 1;
    end
    
    x_pos_sequences(cur_trial,:) = x_filt;
    y_pos_sequences(cur_trial,:) = y_filt;
end

cd ~/Data/bruce/
save drifting_pos_seqs_fin x_pos_sequences y_pos_sequences
