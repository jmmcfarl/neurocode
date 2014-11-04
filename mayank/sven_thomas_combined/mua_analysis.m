clear all
cd C:\WC_Germany\sven_thomas_combined
load combined_dir_nd_dist

%%

for d = 89:length(combined_dir)
    cd(combined_dir{d})
    pwd
    
    if exist('./Sc1.ntt','file')
        
        amp_threshold = 25;
        max_overlap = 0.5;
        [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
        
        load ./sync_times.mat
        synct_d = downsample(synct,8);
        rates = cellfun(@(x) length(x),mua_amps)/range(synct)*1e6;
        %%
        load ./used_data
        clear S
        params.Fs = 2016;
        params.tapers = [2 3];
        params.fpass = [0 40];
        win = 50;
        for i = 2:8
            eval(['[S(i,:),f]=mtspectrumsegc(' sprintf('lf%d',i) ',win,params,1);']);
        end
        uds_freqs = find(f > 0.1 & f < 1);
        peak_uds_pow = max(S(:,uds_freqs),[],2);
        peak_uds_pow(1) = nan;
        
        
        %%
        if exist('./Sc1.ntt','file')
            
            figure
            for i = 1:8
                subplot(8,4,4*(9-i-1)+1)
                plot(avg_waveform(i,:)','linewidth',2)
%             set(gca,'fontsize',20)
                xlim([0 32])
                subplot(8,4,4*(9-i-1)+2)
                hist(mua_amps{i},150)
%             set(gca,'fontsize',20)
                xlim([0 120])
            end
            subplot(8,4,(1:8)*4-1)
            plot([nan rates(2:8)],1:8,'o-','linewidth',2)
%             set(gca,'fontsize',20)
            if ~isnan(hpc_mua(d))
                title(sprintf('MUA used %d',hpc_mua(d)))
            end
            subplot(8,4,(1:8)*4)
            plot(peak_uds_pow,1:8,'o-','linewidth',2)
%             set(gca,'fontsize',20)
            fillPage(gcf,'Papersize',[30 30],'Margins',[0 0 0 0])
        else
            
            plot(peak_uds_pow,1:8,'o-')
            xlabel('UDS Pow','fontsize',16)
            ylabel('Channel','fontsize',16)
            
        end
        temp = find(combined_dir{d} == '\',1,'last');
        cur_name = combined_dir{d}(temp+1:end);
        fname = ['C:\WC_Germany\sven_thomas_combined\MUA\' cur_name];
        print('-djpeg',fname);close all
    end
end
