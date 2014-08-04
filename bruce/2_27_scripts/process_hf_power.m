clear all
close all
cd ~/Data/bruce/2_27_12/M232/

Fs = 40e3; %sample rate
dsf = 40;
Fsd = Fs/dsf;
sm_win = round(Fs*4e-3);

n_ch = 24;
gamma = 1/3; %compressive NL
for blockid = 1:5
    
    fprintf('Block %d of %d\n',blockid,5);
    fname = sprintf('Expt5%dFullV.mat',blockid);
    load(fname);
    
    blklen = FullV.blklen;
    n_blks = length(blklen);
    
    hf_pow = zeros(n_ch,0); %initialize power matrix
    t = [];
    
    cur_point = 1;
    for i = 1:n_blks
        fprintf('Block %d of %d\n',i,n_blks);
        cur_set = cur_point:(cur_point + blklen(i)-1);
        
        temp = zeros(n_ch,blklen(i));
        for ch = 1:n_ch
           temp(ch,:) = double(FullV.V(ch,cur_set)); 
           temp(ch,:) = abs(temp(ch,:)).^gamma; %compress the magnitude of the signal
           temp(ch,:) = smooth(temp(ch,:),sm_win); %boxcar smoothing
        end
        
        hf_pow = [hf_pow; downsample(temp',dsf)];
        t = [t; downsample(FullV.t(cur_set)',dsf)];
        cur_point = cur_point + blklen(i);
    end
    
    sname = sprintf('Expt5%dhfPow',blockid);
    save(sname,'hf_pow','t');
    clear FullV temp hf_pow t
        
end