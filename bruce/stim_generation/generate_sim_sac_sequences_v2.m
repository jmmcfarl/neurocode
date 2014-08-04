close all
clear all

save_dir = '~/Data/bruce/Expt_1_8_13_imfolder';
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

sac_amp = 1.5;
sacs_per_im = 4;

N_nat_trials = 200;
N_noise_trials = 200;
N_flip_trials = 200;
seqs_per_trial = 2;

%% FOR NAT IMAGE TRIALS
nat_dir = '~/James_scripts/data_processing/Images/image_set_A/'
file_range_nat = [1:685];
rr = range(file_range_nat)+1;
file_order_nat = randperm(rr) + file_range_nat(1)-1;

cnt = 1;
trial_offset = 0;
cd(nat_dir)
for cur_trial = 1:N_nat_trials
    fprintf('Nat image trial %d of %d\n',cur_trial,N_nat_trials);
    cur_im = 1;
    for cur_seq = 1:seqs_per_trial
        
        myimgfile = sprintf('%.4d.png',file_order_nat(cnt));
        cnt = cnt + 1;
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        imdataorg = imdata;
        
        %%
        cur_sac_amps = ones(sacs_per_im,1)*sac_amp;
        rand_dirs = rand(sacs_per_im,1)*2*pi-pi;
        X = zeros(sacs_per_im+1,2);
        its = 0;
        for i = 1:sacs_per_im
            X(i+1,1) = X(i,1) + cur_sac_amps(i)*cos(rand_dirs(i));
            X(i+1,2) = X(i,2) + cur_sac_amps(i)*sin(rand_dirs(i));
        end
        X_inds = round(X/Pix2Deg);
        
        %%
        [Ny,Nx] = size(imdata);
        %% cd(save_dir)
        for cur_fix = 1:sacs_per_im
            shift_im = shift_mat2(imdataorg,X_inds(cur_fix,1),X_inds(cur_fix,2),1);
            imdata = shift_im;
            
            %print raw image
            imdatapng = imdata./whiteColor;
            
            im_id = 1000000+1000*(cur_trial+trial_offset) + cur_im;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');
            cur_im = cur_im + 1;
        end
    end
end

%% FOR NOISE TRIALS
noise_dir = save_dir;
file_range_noise = [1:500];
rr = range(file_range_noise)+1;
file_order_noise = randperm(rr) + file_range_noise(1)-1;

trial_offset = N_nat_trials;
cnt = 1;
cd(noise_dir)
for cur_trial = 1:N_noise_trials
    fprintf('Noise image trial %d of %d\n',cur_trial,N_nat_trials);
    cur_im = 1;
    for cur_seq = 1:seqs_per_trial
        cur_id = file_order_noise(cnt);
        myimgfile = sprintf('IM%.7d.png',cur_id);
        cnt = cnt + 1;
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        imdataorg = imdata;
        
        %%
        cur_sac_amps = ones(sacs_per_im,1)*sac_amp;
        rand_dirs = rand(sacs_per_im,1)*2*pi-pi;
        X = zeros(sacs_per_im+1,2);
        its = 0;
        for i = 1:sacs_per_im
            X(i+1,1) = X(i,1) + cur_sac_amps(i)*cos(rand_dirs(i));
            X(i+1,2) = X(i,2) + cur_sac_amps(i)*sin(rand_dirs(i));
        end
        X_inds = round(X/Pix2Deg);
        
        %%
        [Ny,Nx] = size(imdata);
        %% cd(save_dir)
        for cur_fix = 1:sacs_per_im
            shift_im = shift_mat2(imdataorg,X_inds(cur_fix,1),X_inds(cur_fix,2),1);
            imdata = shift_im;
            
            %print raw image
            imdatapng = imdata./whiteColor;
            
            im_id = 1000000+1000*(cur_trial+trial_offset) + cur_im;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');
            cur_im = cur_im + 1;
        end
    end
end

%% FOR NAT FLIP TRIALS
nat_dir = '~/James_scripts/data_processing/Images/image_set_A/'
file_range_nat = [1:685];
rr = range(file_range_nat)+1;
file_order_nat = randperm(rr) + file_range_nat(1)-1;

trial_offset = N_nat_trials + N_noise_trials;
cd(nat_dir)
for cur_trial = 1:N_flip_trials
    fprintf('Flip image trial %d of %d\n',cur_trial,N_nat_trials);
    cur_im = 1;
    for cur_seq = 1:seqs_per_trial*sacs_per_im
        
        cur_im_rand = floor(rand*rr)+file_range_nat(1);
        
        myimgfile = sprintf('%.4d.png',file_range_nat(cur_im_rand));
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
                
        [Ny,Nx] = size(imdata);
        
        %print raw image
        imdatapng = imdata./whiteColor;
        
        im_id = 1000000+1000*(cur_trial+trial_offset) + cur_im;
        cur_fname = sprintf('IM%.7d',im_id);
        imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');
        cur_im = cur_im + 1;
    end
end
