close all
clear all

save_dir = '~/Data/bruce/Expt_1_8_13_imfolder';
blackColor = 0;
whiteColor = 255;
grayColor = (blackColor+whiteColor)/2;
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

cent_x = 0.35;
cent_y = 0.43;
% cent_x = 2;
% cent_y = 2;
% local_rad = 1;
% DD = sqrt((XAX-cent_x).^2 + (YAX-cent_y).^2);
% foreground = DD < local_rad;
local_std = 0.5;
peak_amp = 4;
fore_weight = peak_amp*exp(-((XAX-cent_x).^2 + (YAX-cent_y).^2)/(2*local_std^2));
fore_weight(fore_weight > 1) = 1;
back_weight = 1-fore_weight;

% center pixel before scale transformation
XcenterPix = Nxp/2;
YcenterPix = Nyp/2;

freq_cent = 4.5; %cyc/deg
freq_a = 9;
orange = [-180 180];
dec_const = 2/niqf;

fix_per_seq = 8;
ims_per_fix = 20;
%%
nat_trials = 1:200;
% nat_trials = 1;
%% WITH SIM SAC BACKGND
N_sim_seqs = 3;
seq_offset = 0;
cd(save_dir)
for cur_seq = 1:N_sim_seqs
    cur_im = 1;
    for cur_fix = 1:fix_per_seq/2
        
        cur_id = 1000000+1000*nat_trials(cur_seq) + cur_fix;
        myimgfile = sprintf('IM%.7d.png',cur_id);
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        backgnd1 = imdata;
       
        cur_id = 1000000+1000*nat_trials(cur_seq) + cur_fix + fix_per_seq/2;
        myimgfile = sprintf('IM%.7d.png',cur_id);
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        backgnd2 = imdata;
        
        for ii = 1:ims_per_fix
            
%             [noise_im,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,orange,dec_const);
%             noise_im(Nyp+1:end,:) = [];
%             
%             bounds = prctile(noise_im(:),cbounds);
%             noise_im = (noise_im-bounds(1))/range(bounds);
%             noise_im = noise_im*whiteColor;
%             noise_im(noise_im < blackColor) = blackColor;
%             noise_im(noise_im > whiteColor) = whiteColor;
            
cur_id = 1000*nat_trials(cur_seq) + (cur_fix-1)*ims_per_fix + ii;
myimgfile = sprintf('IM%.7d.png',cur_id);

imdata = imread(myimgfile);
noise_im = double(imdata); % convert to double format

            %             image = backgnd;
            %             image(foreground) = noise_im(foreground);
            image1 = backgnd1.*back_weight + noise_im.*fore_weight;
            image2 = backgnd2.*back_weight + noise_im.*fore_weight;
            
            %print raw image
          imdatapng = image1./whiteColor;            
            im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');
            
            imdatapng = image2./whiteColor;
            im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im + ims_per_fix*fix_per_seq/2;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');
            
            cur_im = cur_im + 1;
        end
    end
end

%% WITH GRAY BACKGND
seq_offset = N_sim_seqs;
N_gray_seqs = 2;
backgnd = ones(Nyp,Nxp)*grayColor;
for cur_seq = 1:N_gray_seqs
    cur_im = 1;
    for ii = 1:ims_per_fix*fix_per_seq/2
        
%         [noise_im,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,orange,dec_const);
%         noise_im(Nyp+1:end,:) = [];
%         
%         bounds = prctile(noise_im(:),cbounds);
%         noise_im = (noise_im-bounds(1))/range(bounds);
%         noise_im = noise_im*whiteColor;
%         noise_im(noise_im < blackColor) = blackColor;
%         noise_im(noise_im > whiteColor) = whiteColor;
 
cur_id = 1000*nat_trials(cur_seq) + (cur_fix-1)*ims_per_fix + ii;
myimgfile = sprintf('IM%.7d.png',cur_id);

imdata = imread(myimgfile);
noise_im = double(imdata); % convert to double format

        %             image = backgnd;
        %             image(foreground) = noise_im(foreground);
        image = backgnd.*back_weight + noise_im.*fore_weight;
        
        %print raw image
        imdatapng = image./whiteColor;
        
        im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im;
        cur_fname = sprintf('IM%.7d',im_id);
        imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');
        
        im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im + ims_per_fix*fix_per_seq/2;
        cur_fname = sprintf('IM%.7d',im_id);
        imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');

        cur_im = cur_im + 1;
    end
end

%% WITH SIM SAC BACKGND (SMALL PATCH)
local_std = 0.25;
peak_amp = 4;
fore_weight = peak_amp*exp(-((XAX-cent_x).^2 + (YAX-cent_y).^2)/(2*local_std^2));
fore_weight(fore_weight > 1) = 1;
back_weight = 1-fore_weight;


N_small_seqs = 2;
seq_offset = N_sim_seqs + N_gray_seqs;
cd(save_dir)
for cur_seq = 1:N_small_seqs
    cur_im = 1;
    for cur_fix = 1:fix_per_seq/2
        
        cur_id = 1000000+1000*nat_trials(cur_seq) + cur_fix;
        myimgfile = sprintf('IM%.7d.png',cur_id);
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        backgnd1 = imdata;
       
        cur_id = 1000000+1000*nat_trials(cur_seq) + cur_fix + fix_per_seq/2;
        myimgfile = sprintf('IM%.7d.png',cur_id);
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        backgnd2 = imdata;
        
        for ii = 1:ims_per_fix
            
%             [noise_im,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,orange,dec_const);
%             noise_im(Nyp+1:end,:) = [];
%             
%             bounds = prctile(noise_im(:),cbounds);
%             noise_im = (noise_im-bounds(1))/range(bounds);
%             noise_im = noise_im*whiteColor;
%             noise_im(noise_im < blackColor) = blackColor;
%             noise_im(noise_im > whiteColor) = whiteColor;

cur_id = 1000*nat_trials(cur_seq) + (cur_fix-1)*ims_per_fix + ii;
myimgfile = sprintf('IM%.7d.png',cur_id);

imdata = imread(myimgfile);
noise_im = double(imdata); % convert to double format

            %             image = backgnd;
            %             image(foreground) = noise_im(foreground);
            image1 = backgnd1.*back_weight + noise_im.*fore_weight;
            image2 = backgnd2.*back_weight + noise_im.*fore_weight;
            
            %print raw image
            imdatapng = image1./whiteColor;            
            im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');
            
            imdatapng = image2./whiteColor;
            im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im + ims_per_fix*fix_per_seq/2;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(save_dir, '/', cur_fname, '.png'),'png');
            
            cur_im = cur_im + 1;
        end
    end
end

%% WITH FLASHED BACKGND
N_flash_seqs = 3;
seq_offset = 7;
nsave_dir = '~/Data/bruce/additional_loc_glob/';
flash_trials = 401:600;
cd(save_dir)
for cur_seq = 1:N_flash_seqs
    cur_im = 1;
    for cur_fix = 1:fix_per_seq/2
        
        cur_id = 1000000+1000*flash_trials(cur_seq) + cur_fix;
        myimgfile = sprintf('IM%.7d.png',cur_id);
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        backgnd1 = imdata;
       
        cur_id = 1000000+1000*flash_trials(cur_seq) + cur_fix + fix_per_seq/2;
        myimgfile = sprintf('IM%.7d.png',cur_id);
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        backgnd2 = imdata;
        
        for ii = 1:ims_per_fix
            
%             [noise_im,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,orange,dec_const);
%             noise_im(Nyp+1:end,:) = [];
%             
%             bounds = prctile(noise_im(:),cbounds);
%             noise_im = (noise_im-bounds(1))/range(bounds);
%             noise_im = noise_im*whiteColor;
%             noise_im(noise_im < blackColor) = blackColor;
%             noise_im(noise_im > whiteColor) = whiteColor;
            
cur_id = 1000*nat_trials(cur_seq) + (cur_fix-1)*ims_per_fix + ii;
myimgfile = sprintf('IM%.7d.png',cur_id);

imdata = imread(myimgfile);
noise_im = double(imdata); % convert to double format

            %             image = backgnd;
            %             image(foreground) = noise_im(foreground);
            image1 = backgnd1.*back_weight + noise_im.*fore_weight;
            image2 = backgnd2.*back_weight + noise_im.*fore_weight;
            
            %print raw image
          imdatapng = image1./whiteColor;            
            im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(nsave_dir, '/', cur_fname, '.png'),'png');
            
            imdatapng = image2./whiteColor;
            im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im + ims_per_fix*fix_per_seq/2;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(nsave_dir, '/', cur_fname, '.png'),'png');
            
            cur_im = cur_im + 1;
        end
    end
end

%% WITH FLASHED BACKGND AND GRATING FORE
N_flash_seqs = 3;
seq_offset = 10;
nsave_dir = '~/Data/bruce/additional_loc_glob/';
flash_trials = 401:600;
cd(save_dir)
for cur_seq = 1:N_flash_seqs
    cur_im = 1;
    for cur_fix = 1:fix_per_seq/2
        
        cur_id = 1000000+1000*flash_trials(cur_seq) + cur_fix;
        myimgfile = sprintf('IM%.7d.png',cur_id);
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        backgnd1 = imdata;
       
        cur_id = 1000000+1000*flash_trials(cur_seq) + cur_fix + fix_per_seq/2;
        myimgfile = sprintf('IM%.7d.png',cur_id);
        
        imdata = imread(myimgfile);
        imdata = double(imdata); % convert to double format
        
        backgnd2 = imdata;
        
        for ii = 1:ims_per_fix
            
            %             [noise_im,w2,Fx,Fy] = create_bandexpfilt_noise(siz,freq_cent/niqf,freq_a/niqf,orange,dec_const);
            %             noise_im(Nyp+1:end,:) = [];
            %
            %             bounds = prctile(noise_im(:),cbounds);
            %             noise_im = (noise_im-bounds(1))/range(bounds);
            %             noise_im = noise_im*whiteColor;
            %             noise_im(noise_im < blackColor) = blackColor;
            %             noise_im(noise_im > whiteColor) = whiteColor;
            
            cur_id = 1000*(nat_trials(cur_seq)+3) + (cur_fix-1)*ims_per_fix + ii;
            myimgfile = sprintf('IM%.7d.png',cur_id);
            
            imdata = imread(myimgfile);
            noise_im = double(imdata); % convert to double format
            noise_im(1:128,:) = []; noise_im(end-127:end,:) = [];
            
            %             image = backgnd;
            %             image(foreground) = noise_im(foreground);
            image1 = backgnd1.*back_weight + noise_im.*fore_weight;
            image2 = backgnd2.*back_weight + noise_im.*fore_weight;
            
            %print raw image
          imdatapng = image1./whiteColor;            
            im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(nsave_dir, '/', cur_fname, '.png'),'png');
            
            imdatapng = image2./whiteColor;
            im_id = 2000000+1000*(cur_seq+seq_offset) + cur_im + ims_per_fix*fix_per_seq/2;
            cur_fname = sprintf('IM%.7d',im_id);
            imwrite(imdatapng,strcat(nsave_dir, '/', cur_fname, '.png'),'png');
            
            cur_im = cur_im + 1;
        end
    end
end

