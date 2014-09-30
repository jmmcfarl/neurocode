clear all
close all

% Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
% ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];

for ee = 1:length(Expt_list)
    for oo = 1:2
        Expt_name = Expt_list{ee};
        bar_ori = ori_list(ee,oo);
        if ~isnan(bar_ori)
            fit_unCor = 0;
            include_bursts = 0;
            
            fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';
            base_tname = 'sac_trig_avg_data';
            base_sname = 'sacStimProcFinR2';
            base_yname = 'sacTypeDep';
            
            if include_bursts
                base_tname = strcat(base_tname,'_withbursts');
                base_sname = strcat(base_sname,'_withbursts');
                base_yname = strcat(base_yname,'_withbursts');
            end
            
            if Expt_name(1) == 'M'
                rec_type = 'LP';
                good_coils = [1 1]; %which coils are usable
                use_coils = [1 1]; %[L R] Use info from coils?
                n_probes = 24;
            elseif Expt_name(1) == 'G'
                good_coils = [1 0]; %which coils are usable
                use_coils = [0 0]; %[L R] Use info from coils?
                n_probes = 96;
            end
            
            et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
            et_mod_data_name = 'full_eyetrack_initmods_Rinit';
            et_anal_name = 'full_eyetrack_Rinit';
            
            %if using coil info
            if any(use_coils > 0)
                et_anal_name = [et_anal_name '_Cprior'];
            end
            
            et_mod_data_name = [et_dir et_mod_data_name sprintf('_ori%d',bar_ori)];
            et_anal_name = [et_dir et_anal_name sprintf('_ori%d',bar_ori)];
            load(et_anal_name,'et_params');
            
            Expt_num = str2num(Expt_name(2:end));
            sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
            
            %load trig avg data
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',bar_ori));
            load(tname);
            %load stimulus proc data
            sname = strcat(sac_dir,base_sname,sprintf('_ori%d',bar_ori));
            if fit_unCor
                sname = strcat(sname,'_unCor');
            end
            load(sname);
            %load type-dep data (STILL NEED TO ADD THIS)
            yname = strcat(sac_dir,base_yname,sprintf('_ori%d',bar_ori));
            if fit_unCor
                yname = strcat(yname,'_unCor');
            end
            load(yname);
            
            su_range = (n_probes+1):length(sacStimProc);
            
            for cc = su_range
                su_ind = find(su_range == cc);
                clear SU_data
                cur_sacStimProc = sacStimProc(su_range(su_ind));
                trig_avg = sua_data(su_ind);
                type_dep = sacTypeDep(su_range(su_ind));
                
                if ~isempty(cur_sacStimProc.gsac_avg_rate)
                    cur_GQM = cur_sacStimProc.ModData.rectGQM;
                    flen = cur_GQM.stim_params(1).stim_dims(1);
                    sp_dx = et_params.sp_dx;
                    use_nPix = et_params.use_nPix;
                    use_nPix_us = use_nPix*et_params.spatial_usfac;
                    
                    tlags = trig_avg_params.lags*trig_avg_params.dt;
                    cid = sprintf('E%d_C%d_',Expt_num,cc);
                    
                    %% STA figs
                    sac_xr = [-0.1 0.3];
                    
                    close all
                    stim_dims = cur_sacStimProc.ModData.rectGQM.stim_params(1).stim_dims;
                    pix_ax = (1:stim_dims(2))*sp_dx;
                    pix_ax = pix_ax - mean(pix_ax);
                    lag_ax = ((1:stim_dims(1))*dt - dt/2)*1e3;
                    
                    ov_sta = cur_sacStimProc.ov_phaseDep_sta;
                    [~,sta_peakloc] = max(sum(ov_sta.^2,2));
                    % sta_peakloc = 7;
                    cond_STA = cur_sacStimProc.gsac_phaseDep_sta;
                    cond_STA = squeeze(cond_STA(:,sta_peakloc,:));
                    
                    space_sm = 0.75;
                    time_sm = 0.75;
                    
                    %smooth conditional STA in space
                    if space_sm > 0
                        for iii = 1:size(cond_STA,2)
                            cond_STA(:,iii) = jmm_smooth_1d_cor(cond_STA(:,iii),space_sm);
                        end
                    end
                    
                    %smooth conditional STA in time
                    if time_sm > 0
                        for iii = 1:size(cond_STA,1)
                            cond_STA(iii,:) = jmm_smooth_1d_cor(cond_STA(iii,:),time_sm);
                        end
                    end
                    
%                     xl = [-0.3 0.3];
                    
                    f1 = figure();
                    subplot(2,2,[3 4])
                    imagesc(slags*dt,pix_ax,cond_STA');
                    cam = max(abs(cond_STA(:)));
                    caxis([-cam cam]*0.9);
                    yl = ylim();
                    xlim(sac_xr)
                    xlabel('Time (s)');
                    ylabel('Rel Position (deg)');
%                     ylim(xl);
                    
                    %             f2 = figure();
                    subplot(2,2,1);
                    imagesc(pix_ax,lag_ax,ov_sta);
                    cam = max(abs(ov_sta(:)));
                    caxis([-cam cam]);
                    cid = sprintf('E%d_C%d_',Expt_num,cc);
                    set(gca,'ydir','normal');
%                     xlim(xl);
%                     line(xl,lag_ax([sta_peakloc sta_peakloc]),'color','k');
                    
                    fig_width = 6;
                    rel_height = 0.75;
                    fname = [fig_dir 'STA_examps/' cid sprintf('ori%d_',bar_ori) 'staexamp.pdf'];
                    exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
                    close(f1);
                    
                end
            end
        end
    end
end