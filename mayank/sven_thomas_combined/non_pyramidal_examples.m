clear all
close all

cur_dir_letter = 'C';
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\
% load ./combined_dir.mat
load ./combined_dir_nd.mat

%%
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.01;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

%% EXAMPLE LECL3 non-pyramidal
d = 63; 

% xl = [22 38.5];
xl = [779.5 798];
cur_dir = combined_dir{d};
    cd(cur_dir)
    pwd
    conv_fact = 2.4426e-5*100/50 %based on rec notes
    
    load ./used_data lf8 lf7 wcv
    if ctx_lfp(d) == 7
        lf8 = lf7;
    elseif ctx_lfp(d) == 6
        lf8 = lf6;
    end
    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv = wcv*conv_fact*1e3; %in mV
    scale = std(wcv);
    wcv = zscore(wcv);
    t_axis2 = (1:length(wcv))/2016;
    twenty = 20/scale;
    
    load ./aligned_heka.mat
    
    plot(t_axis2-xl(1),wcv,'b')
    hold on
    plot(t_axis-xl(1),lf8_lf,'k')
    xlim(xl-xl(1))
    set(gca,'fontsize',14,'fontname','arial')
    xlabel('Time (s)','fontsize',16)
    ylabel('Amplitude (z)','fontsize',16)
    line([3 3],[-1 (-1+twenty)],'color','k')
    
    
    %% EXAMPLE MECL3 non-pyramidal
    xl = [26.5 44];
    d = 57; %57 [600 - 625]
cur_dir = combined_dir{d};
    cd(cur_dir)
    pwd
    conv_fact = 2.4426e-5*100/90 %based on rec notes
    
    load ./used_data lf8 lf7 wcv
    if ctx_lfp(d) == 7
        lf8 = lf7;
    elseif ctx_lfp(d) == 6
        lf8 = lf6;
    end
    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv = wcv*conv_fact*1e3; %in mV
    scale = std(wcv);
    wcv = zscore(wcv);
    t_axis2 = (1:length(wcv))/2016;
    twenty = 20/scale;
    
%     load ./aligned_heka.mat
    
    plot(t_axis2-xl(1),wcv,'b')
    hold on
    plot(t_axis-xl(1),lf8_lf,'k')
    xlim(xl-xl(1))
    set(gca,'fontsize',14,'fontname','arial')
    xlabel('Time (s)','fontsize',16)
    ylabel('Amplitude (z)','fontsize',16)
    line([3 3],[-1 (-1+twenty)],'color','k')
    
    %% EXAMPLE MECL3 Interneuron
d = 62;
xl = [116 132.5];
cur_dir = combined_dir{d};
    cd(cur_dir)
    pwd
    conv_fact = 2.4426e-5*100/80 %based on rec notes
    
    load ./used_data lf8 lf7 wcv
    if ctx_lfp(d) == 7
        lf8 = lf7;
    elseif ctx_lfp(d) == 6
        lf8 = lf6;
    end
    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv = wcv*conv_fact*1e3; %in mV
    scale = std(wcv);
    wcv = zscore(wcv);
    t_axis2 = (1:length(wcv))/2016;
    twenty = 20/scale;
    
%     load ./aligned_heka.mat
    
    plot(t_axis2-xl(1),wcv,'b')
    hold on
    plot(t_axis-xl(1),lf8_lf,'k')
    xlim(xl-xl(1))
    set(gca,'fontsize',14,'fontname','arial')
    xlabel('Time (s)','fontsize',16)
    ylabel('Amplitude (z)','fontsize',16)
    line([3 3],[-1 (-1+twenty)],'color','k')
