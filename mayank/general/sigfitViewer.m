%sigmoid fit viewer

%for wcv
sigFitWin = [-800:800];

for i = 5:length(wcv_up_start_id)
    begID = wcv_up_start_id(i) - round(5/dt);
    endID = wcv_up_start_id(i) + round(5/dt);
    taxis = [begID:endID]*dt;
    offset = min(wcv(wcv_up_start_id(i)+sigFitWin));
    plot(taxis,wcv(begID:endID),'linewidth',2)
    hold on
    plot((sigFitWin+wcv_up_start_id(i))*dt,my_sigmoid(net_betafit_up(i,:),sigFitWin)+offset,'r','linewidth',2)
    
    offset = max(wcv(wcv_up_end_id(i)+sigFitWin));
    plot((sigFitWin+wcv_up_end_id(i))*dt,-my_sigmoid(net_betafit_down(i,:),sigFitWin)+offset,'g','linewidth',2)
    
    pause
    clf
end