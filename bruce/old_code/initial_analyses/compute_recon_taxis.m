clear all
close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat

block_id = 1;
block_times = Blocks{block_id}.blocktimes;
stim_times = Blocks{block_id}.stimtime;

%%
dt = 0.05;

cd stimrecon/
recon_t = [];
stim_num = [];
for i = 1:length(stim_times)
    i
    load(['BLOCK' num2str(block_id) 'IMAGE' num2str(i) '.mat'], 'STIMrec');
    
    n_recon_samps(i) = size(STIMrec,1);
    recon_t = [recon_t stim_times(i):dt:(stim_times(i)+n_recon_samps(i)*dt)];
    stim_num = [stim_num i*ones(1,n_recon_samps(i))];
end

cd ~/Data/bruce/2_27_12/stimrecon/
save stimrecon_t dt stim_num recon_t n_recon_samps