% close all

dir_prefix = '~';
% dir_prefix = '/media/NTlab_data1';
Expt_name = 'G088';
data_dir = ['~/Data/bruce/' Expt_name];
raw_data_dir = [dir_prefix '/Data/bruce/' Expt_name];

% save_dir = data_dir;
save_dir = ['~/Data/bruce/' Expt_name];

cd(data_dir);
load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

force_recompute = 1;
%%
use_lfps = 1:96;
Fs = 1/3.333307000208032e-05;
dsf = 75;
Fsd = Fs/dsf;
niqf = Fsd/2;
hcf = 0.8*niqf;
filt_ord = 4;
[filt_b,filt_a] = butter(filt_ord,hcf/(Fs/2),'low');

lfp_params.Fs = Fs;
lfp_params.Fsd = Fsd;
lfp_params.dsf = dsf;
lfp_params.hcf = hcf;
lfp_params.filt_type = 'butter';
lfp_params.filt_ord = filt_ord;

cd(raw_data_dir);
for ee = 1:length(Expts)
    
    d = dir(sprintf('Expt%d.p*FullV.mat',ee));
    fname = [save_dir sprintf('/Expt%d_LFP.mat',ee)];
    lfp_exist = exist(fname,'file');
    if force_recompute
        lfp_exist = 0;
    end
    if ~isempty(d) & ~lfp_exist
        fprintf('Processing Expt %d of %d\n',ee,length(Expts));
        
        filename = sprintf('Expt%dFullVmean.mat',ee);
        load(filename);
        
        lfp_mat = [];
        lfp_int2V = nan(length(use_lfps),1);
        for ll = 1:length(use_lfps)
            fprintf('Electrode %d of %d\n',ll,length(use_lfps));
            filename = sprintf('Expt%d.p%dFullV.mat',ee,use_lfps(ll));
            load(filename);
            V = double(FullV.V);
            lfp_int2V(ll) = FullV.intscale(1)/FullV.intscale(2);
            
            V = V + FullV.sumscale*sumv;
            V = V*lfp_int2V(ll);
            nparts = length(FullV.blklen);
            cur_lfp = [];
            %splice together multiple blocks
            cur_pt = 1;
            for pp = 1:nparts
                cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
                cur_range(cur_range > length(V)) = [];
                cur_V = filtfilt(filt_b,filt_a,V(cur_range));
                cur_V = downsample(cur_V,dsf);
                cur_lfp = [cur_lfp cur_V];
                cur_pt = cur_pt + FullV.blklen(pp);
            end
            cur_lfp = int16(cur_lfp'/lfp_int2V(ll));
            lfp_mat = cat(2,lfp_mat,cur_lfp);
        end
        
        lfp_t_ax = [];
        for pp = 1:nparts
            cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
            lfp_t_ax = [lfp_t_ax downsample(cur_t_ax,dsf)];
        end
        if length(lfp_t_ax) > size(lfp_mat,1)
            fprintf('%d extra time points, removing\n',length(lfp_t_ax) - size(lfp_mat,1));
            lfp_t_ax(size(lfp_mat,1)+1:end) = [];
        end
        
        fprintf('LFP len: %d sec\n',range(lfp_t_ax));
        
        save(fname,'lfp_t_ax','lfp_mat','lfp_params','lfp_int2V');
    elseif lfp_exist
        fprintf('LFP file exists, skipping\n');        
    else
        fprintf('No files found for expt %d, skipping\n',ee);
    end
end
