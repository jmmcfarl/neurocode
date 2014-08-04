clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\overall_info_file_data
load C:\WC_Germany\overall_calcs\eyeball_theta_times

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.fpass = [0 10];
params.err = [2 .05];
movingwin = [30 30];

niqf = 2016/2;
hcf1 = 100/niqf;
[b,a] = butter(2,hcf1,'low');

for d = 1:33

    cd(over_dir{d});
    disp(['session ' num2str(d)]);

    load used_data lf8 lf3 wcv_minus_spike

    down_w = downsample(wcv_minus_spike/mp_gain(d),dsf);
    down_8 = downsample(lf8/lf8_gain(d),dsf);
    down_3 = downsample(lf3/lf3_gain(d),dsf);

    if ~isempty(theta_times{d}) && ~isempty(nt_times{d})

        tot_theta = sum(theta_times{d}(:,2) - theta_times{d}(:,1));
        tot_nt = sum(nt_times{d}(:,2) - nt_times{d}(:,1));
        theta_perc(d) = tot_theta/(tot_theta+tot_nt);

        [Sw_theta(d,:), f, Swerr_theta(d,:,:)]= mtspectrumc_unequal_length_trials...
            (down_w, movingwin, params,theta_times{d}*Fsd+1);
        [Sw_nt(d,:), f, Swerr_nt(d,:,:)]= mtspectrumc_unequal_length_trials...
            (down_w, movingwin, params,nt_times{d}*Fsd+1);

        [S8_theta(d,:), f, S8err_theta(d,:,:)]= mtspectrumc_unequal_length_trials...
            (down_8, movingwin, params,theta_times{d}*Fsd+1);
        [S8_nt(d,:), f, S8err_nt(d,:,:)]= mtspectrumc_unequal_length_trials...
            (down_8, movingwin, params,nt_times{d}*Fsd+1);

        [S3_theta(d,:), f, S3err_theta(d,:,:)]= mtspectrumc_unequal_length_trials...
            (down_3, movingwin, params,theta_times{d}*Fsd+1);
        [S3_nt(d,:), f, S3err_nt(d,:,:)]= mtspectrumc_unequal_length_trials...
            (down_3, movingwin, params,nt_times{d}*Fsd+1);


        [Cw8_theta(d,:),Phiw8_theta(d,:),Smn,Smm,f,ConfCw8_theta(d),PhiStd,Cerr] = coherencyc_unequal_length_trials...
            ([down_w down_8], movingwin, params, theta_times{d}*Fsd+1);
        [Cw8_nt(d,:),Phiw8_nt(d,:),Smn,Smm,f,ConfCw8_nt(d),PhiStd,Cerr] = coherencyc_unequal_length_trials...
            ([down_w down_8], movingwin, params, nt_times{d}*Fsd+1);

        [Cw3_theta(d,:),Phiw3_theta(d,:),Smn,Smm,f,ConfCw3_theta(d),PhiStd,Cerr] = coherencyc_unequal_length_trials...
            ([down_w down_3], movingwin, params, theta_times{d}*Fsd+1);
        [Cw3_nt(d,:),Phiw3_nt(d,:),Smn,Smm,f,ConfCw3_nt(d),PhiStd,Cerr] = coherencyc_unequal_length_trials...
            ([down_w down_3], movingwin, params, nt_times{d}*Fsd+1);

        [C83_theta(d,:),Phi83_theta(d,:),Smn,Smm,f,ConfC83_theta(d),PhiStd,Cerr] = coherencyc_unequal_length_trials...
            ([down_8 down_3], movingwin, params, theta_times{d}*Fsd+1);
        [C83_nt(d,:),Phi83_nt(d,:),Smn,Smm,f,ConfC83_nt(d),PhiStd,Cerr] = coherencyc_unequal_length_trials...
            ([down_8 down_3], movingwin, params, nt_times{d}*Fsd+1);


        plot(f,10*log10(Sw_theta(d,:)))
        hold on
        plot(f,10*log10(Sw_nt(d,:)),'r')
        title(['Percent Theta: ' num2str(theta_perc(d))]);
        tname = ['C:\WC_Germany\overall_calcs\state_dep_spec\MP_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close all

        plot(f,10*log10(S8_theta(d,:)))
        hold on
        plot(f,10*log10(S8_nt(d,:)),'r')
        title(['Percent Theta: ' num2str(theta_perc(d))]);
        tname = ['C:\WC_Germany\overall_calcs\state_dep_spec\LF8_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close all

        plot(f,10*log10(S3_theta(d,:)))
        hold on
        plot(f,10*log10(S3_nt(d,:)),'r')
        title(['Percent Theta: ' num2str(theta_perc(d))]);
        tname = ['C:\WC_Germany\overall_calcs\state_dep_spec\LF3_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close all

        plot(f,Cw8_theta(d,:))
        hold on
        plot(f,Cw8_nt(d,:),'r')
        title(['Percent Theta: ' num2str(theta_perc(d))]);
        line([0 max(f)],[ConfCw8_theta(d) ConfCw8_theta(d)])
        line([0 max(f)],[ConfCw8_nt(d) ConfCw8_nt(d)],'Color','r')
        tname = ['C:\WC_Germany\overall_calcs\state_dep_spec\Cw8_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close all

        plot(f,Cw3_theta(d,:))
        hold on
        plot(f,Cw3_nt(d,:),'r')
        title(['Percent Theta: ' num2str(theta_perc(d))]);
        line([0 max(f)],[ConfCw3_theta(d) ConfCw3_theta(d)])
        line([0 max(f)],[ConfCw3_nt(d) ConfCw3_nt(d)],'Color','r')
        tname = ['C:\WC_Germany\overall_calcs\state_dep_spec\Cw3_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close all

        plot(f,C83_theta(d,:))
        hold on
        plot(f,C83_nt(d,:),'r')
        title(['Percent Theta: ' num2str(theta_perc(d))]);
        line([0 max(f)],[ConfC83_theta(d) ConfC83_theta(d)])
        line([0 max(f)],[ConfC83_nt(d) ConfC83_nt(d)],'Color','r')
        tname = ['C:\WC_Germany\overall_calcs\state_dep_spec\C83_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close all

    end

end


