% eye calibration
clear all
datafolder = './Data/ExptB/';
savetofolder = './Data/ExptBproc/';
% addpath('./Data/ExptA/'); load EyeCal17Feb2012
%  eyecal = EyeCalD.Data;
addpath(datafolder); load EyeCal20feb2012
eyecal = EyeCalData.Data;

eyecal(5,:) = -eyecal(5,:); % flip x-location


FixXu = unique(eyecal(5,:));
FixYu = unique(eyecal(6,:));

FixXu = FixXu(~isnan(FixXu));
FixYu = FixYu(~isnan(FixYu));
FixXu = FixXu(abs(FixXu)<7);
FixYu = FixYu(abs(FixYu)<7);
EyeCluster = [];
EyeClusterStd = [];
FixCenter = [];
EyeLdatas = [];
EyeRdatas = [];

for i=1:length(FixXu)
    for j=1:length(FixYu)
        inds = find(eyecal(5,:)==FixXu(i) & eyecal(6,:)==FixYu(j));
        if isempty(inds)
            continue;
        end
        
        eyecalij = eyecal(1:4, inds);
        centerij = median(eyecalij');
        Ndata = size(eyecalij,2);
        dists = abs(eyecalij-centerij'*ones(1, Ndata));
        
        goodpt = find(any(dists<.5)==1);
        eyecalijgood = eyecalij(:, goodpt);
        
        Ngoodpt = length(goodpt);
        boundLOW = prctile(eyecalijgood', 20);
        boundHIGH = prctile(eyecalijgood', 80);
        
        InRange = eyecalijgood>boundLOW'*ones(1,Ngoodpt) & eyecalijgood<boundHIGH'*ones(1,Ngoodpt);
        EyeLpt = find(and(InRange(1, :), InRange(3, :))>0);
        EyeRpt = find(and(InRange(2, :), InRange(4, :))>0);
        
        CenterEst = zeros(1,4);
        NoiseEst = zeros(1,4);
        for k=1:4
            data = eyecalijgood(k,:);
%             boundLOW = prctile(data, 20);
%             boundHIGH = prctile(data, 80);
%             data = data(data>boundLOW & data<boundHIGH);
            if k==1 || k==3
                pts = EyeLpt;
            elseif k==2 || k==4
                pts = EyeRpt;
            end
            CenterEst(k) = mean(data(pts));
            NoiseEst(k) = std(data(pts));            
        end
        EyeLdatas = [EyeLdatas; [eyecalijgood([1 3], EyeLpt)' ones(size(EyeLpt))'*[FixXu(i) FixYu(j)] ] ]; 
        EyeRdatas = [EyeRdatas; [eyecalijgood([2 4], EyeRpt)' ones(size(EyeRpt))'*[FixXu(i) FixYu(j)] ]]; 
        FixCenter = [FixCenter; [FixXu(i) FixYu(j)]];
        
        EyeCluster = [EyeCluster; CenterEst];
        EyeClusterStd = [EyeClusterStd NoiseEst];
    end
end

gains = zeros(4, 2);
offsets = zeros(4, 1);
highergain = zeros(4, 5);
% funObj = @(p, x) x*p(1:2)'+p(3);
costObj = @(p, x, y) sum((x*p(1:2)+p(3) - y).^2);
p0 = [1 1 0]';
% costObj = @(p, x, y) sum((x*p(1:2)+p(3) + p(4)*x(:,1).^2+ p(5)*x(:,2).^2+ p(6)*x(:,1).*x(:,2) - y).^2);
% p0 = [1 1 0 0 0 0]';
opts = optimset('GradObj','off','Display','off','MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-7);
for eyei=1:2
    for hv=1:2
        EyeMeasure = EyeCluster(:,[eyei eyei+2]);
        
%         beta = nlinfit(EyeMeasure,FixCenter(:,hv),funObj,[-1 0 0]);
        if eyei==1
            beta = fminsearch(@(p)costObj(p, EyeLdatas(:,1:2),EyeLdatas(:,2+hv)),p0,opts);
        elseif eyei==2
            beta = fminsearch(@(p)costObj(p, EyeRdatas(:,1:2),EyeRdatas(:,2+hv)),p0,opts);
        end
        gains(2*(eyei-1)+hv,:) = beta(1:2);
        offsets(2*(eyei-1)+hv,:) = beta(3);
        if length(beta)>3
            highergain(2*(eyei-1)+hv,1:(length(beta)-3)) = beta(4:end);
        end
%                 figure,plot(funObj(beta, EyeMeasure), FixCenter(:,hv),'r.')
    end
end

% after calibration
EyeLdatascal = CalibrateEyeSignal(EyeLdatas(:,1:2), gains(1:2,:),offsets(1:2), highergain(1:2,:));
EyeRdatascal = CalibrateEyeSignal(EyeRdatas(:,1:2), gains(3:4,:),offsets(3:4), highergain(3:4,:));

% estimate error
rMSEs = zeros(4,length(FixXu), length(FixYu));
medianDIST = zeros(4,length(FixXu), length(FixYu));
for i=1:length(FixXu)
    for j=1:length(FixYu)
        % left
        Linds = find(EyeLdatas(:,3)==FixXu(i) & EyeLdatas(:,4)==FixYu(j));
        rMSEs(1,i,j) = sqrt(median((EyeLdatascal(Linds,1) - FixXu(i)).^2));
        rMSEs(2,i,j) = sqrt(median((EyeLdatascal(Linds,2) - FixYu(j)).^2));
        medianDIST(1,i,j) = abs(median(EyeLdatascal(Linds,1)) - FixXu(i));
        medianDIST(2,i,j) = abs(median(EyeLdatascal(Linds,2)) - FixYu(j));
        % right
        Rinds = find(EyeRdatas(:,3)==FixXu(i) & EyeRdatas(:,4)==FixYu(j));
        rMSEs(3,i,j) = sqrt(median((EyeRdatascal(Rinds,1) - FixXu(i)).^2));
        rMSEs(4,i,j) = sqrt(median((EyeRdatascal(Rinds,2) - FixYu(j)).^2));  
        medianDIST(3,i,j) = abs(median(EyeRdatascal(Rinds,1)) - FixXu(i));
        medianDIST(4,i,j) = abs(median(EyeRdatascal(Rinds,2)) - FixYu(j));        
    end
end
figure,
subplot(2,2,1);imagesc(squeeze(rMSEs(1,:,:)));colorbar;title('L Eye-H');
subplot(2,2,2);imagesc(squeeze(rMSEs(2,:,:)));colorbar;title('L Eye-V');
subplot(2,2,3);imagesc(squeeze(rMSEs(3,:,:)));colorbar;title('R Eye-H');
subplot(2,2,4);imagesc(squeeze(rMSEs(4,:,:)));colorbar;title('R Eye-V');
figure,
subplot(2,2,1);imagesc(squeeze(medianDIST(1,:,:)));colorbar;title('L Eye-H');
subplot(2,2,2);imagesc(squeeze(medianDIST(2,:,:)));colorbar;title('L Eye-V');
subplot(2,2,3);imagesc(squeeze(medianDIST(3,:,:)));colorbar;title('R Eye-H');
subplot(2,2,4);imagesc(squeeze(medianDIST(4,:,:)));colorbar;title('R Eye-V');
mean(mean(mean(rMSEs(1:2,:,:))))
mean(mean(mean(rMSEs(3:4,:,:))))
mean(mean(mean(medianDIST(1:2,:,:))))
mean(mean(mean(medianDIST(3:4,:,:))))
% save data
EyePara = struct('gain', gains, 'offset', offsets);
EyePara.label = {'LeyeH','LeyeV','ReyeH','ReyeV'};
EyePara.gain
EyePara.offset

save([savetofolder 'EyeCalParas'],'EyePara');

% show data
axissize = [FixXu(1) FixXu(end) FixYu(1) FixYu(end)]*1.2;
h=figure(1);clf
subplot(2,2,1),plot(EyeLdatas(:,1), EyeLdatas(:,2),'.');hold on;
plot(FixCenter(:,1), FixCenter(:,2),'r.','MarkerSize',12);axis square
title('Before Calibration-L');axis(axissize);

subplot(2,2,2),plot(EyeRdatas(:,1), EyeRdatas(:,2),'.');hold on;
plot(FixCenter(:,1), FixCenter(:,2),'r.','MarkerSize',12);axis square
title('Before Calibration-R');axis(axissize);
subplot(2,2,3),plot(EyeLdatascal(:,1), EyeLdatascal(:,2),'.');hold on;
plot(FixCenter(:,1), FixCenter(:,2),'r.','MarkerSize',12);axis square
title('After Calibration-L');axis(axissize);

subplot(2,2,4),plot(EyeRdatascal(:,1), EyeRdatascal(:,2),'.');hold on;
plot(FixCenter(:,1), FixCenter(:,2),'r.','MarkerSize',12);axis square
title('After Calibration-R');axis(axissize);

print(h,'-dpdf',[savetofolder 'EyeCalibration.pdf']);