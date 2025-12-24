clear; close all; clc;
tic
fsize=16;
%% Show method of putting an exponential lag onto data
dz=.25;
depth=20; % As d and thus t go to infinity, lag goes to tau.
z=depth:-dz:dz;

t=1:length(z); % seconds
fakeT = 8+4*t/max(t);
fakeT(9.5<fakeT & fakeT<=10.5)=10;
tauT=2; % response time 63% (1/e) is 2 sec.
lagtT=lagTime(t,tauT);
lagfakeT=interp1(lagtT,fakeT,t);

fakeO = 0+300*t/max(t);
fakeO(100<fakeO & fakeO<=150)=125;
tauO=25; % response time 63% (1/e) is 25 sec.
lagtO=lagTime(t,tauO);
lagfakeO=interp1(lagtO,fakeO,t);

figure(100);
subplot(221)
plot(fakeT,z,lagfakeT,z);
set(gca,'ydir','rev');
ylabel('depth (m)'); xlabel('Temperature (Celsius)'); legend({'T','Lagged T'},'Location','southeast');
title(['Temperature Tau is ',int2str(tauT),' seconds']); grid on;
subplot(222)
plot(t,t,t,lagtT); xlabel('time'); ylabel('time'); axis equal;
legend({'t','Lagged t'},'Location','southeast');
grid on; axis(gca,[0 max(t) 0 max(t)]);

subplot(223)
plot(fakeO,z,lagfakeO,z);
set(gca,'ydir','rev');
ylabel('depth (m)'); xlabel('Temperature (Celsius)'); legend({'DO','Lagged DO'},'Location','southeast');
title(['Oxygen Tau is ',int2str(tauO),' seconds']); grid on;
subplot(224)
plot(t,t,t,lagtO); xlabel('time'); ylabel('time'); axis equal;
legend({'t','Lagged t'},'Location','southeast');
grid on; axis(gca,[0 max(t) 0 max(t)]);


%% test lags with data
folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder); load CE01ISSPdosta.mat; CTD=load('CE01ISSP.mat');
% [ce.Properties.VariableNames',ce.Properties.VariableUnits']

% do for CE01 
% DONE. Brandy's method is consistent with mine. "CSPP DO data are averaged from 2-7 meters depth".
% DONE. check salinity spikes. See if I can correct by shifting conductivity and recalculating. 
%       No. Spikes are not from a conductivity lag. Bubbles.
% DONE. Optode T is lagged. See if the method fixes it. Yes.
% DONE. Show with improved r2.
% Find max of time lag covariance function. It's closer to 3 sec than 2 sec for optode temperature because of flow distortion.

% plot T-DO. Note problem that thermocline is at the wrong depth
% fix DO. Show fixed DO puts the oxycline at the right depth on poster

for i=23 %1:max(ce.deployment)
    dex = find(ce.deployment==i);
    if isempty(dex)
        nsite=1;
        disp(['nsite ',int2str(nsite),' no data deployment ',int2str(i)]);
    else
        depth = ce.depth(dex)-.50; % vertical offset from length of CTD
        Time = ce.Time(dex);
        %sea_water_electrical_conductivity = ce.sea_water_electrical_conductivity(dex);
        sea_water_practical_salinity = ce.sea_water_practical_salinity(dex);
        sea_water_temperature = ce.sea_water_temperature(dex);
        optode_temperature = ce.optode_temperature(dex);
        dissolved_oxygen = ce.dissolved_oxygen(dex); % umol kg-1
        optode_oxygen = ce.dosta_abcdjm_cspp_tc_oxygen(dex); % umol L-1
        % estimated_oxygen_concentration = double(ce.estimated_oxygen_concentration(dex)); % umol L-1 SAME AS dosta_abcdjm_cspp_tc_oxygen
        profile = ce.profile(dex);
        lagT=optode_oxygen*NaN;
        lagO=lagT;

        % make a mean DO - T relationships
        figure(102);
        dex=find(depth<7);
        plot(optode_temperature(dex),optode_oxygen(dex),'mx','LineWidth',2);
        xlabel('temperature (Celsius)'); ylabel('dissolved oxygen (uM)'); hold on; grid on;
        dex=find(depth>15);
        plot(optode_temperature(dex),optode_oxygen(dex),'gx','LineWidth',2);
        % get MFN and NSIF values
        ind = find( min(Time) < bottom.Time & bottom.Time < max(Time) );
        rnd = find( min(Time) <  riser.Time &  riser.Time < max(Time) );
        plot(bottom.sea_water_temperature(ind),bottom.ctd_tc_oxygen(ind),'.r');
        plot( riser.sea_water_temperature(rnd), riser.ctd_tc_oxygen(rnd),'.b');

        % plot all lines with no values over sat100
        sat100=600; % https://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
        x=optode_temperature;               y=optode_oxygen; 
        dex=find(y<sat100 & depth<7); x=x(dex); y=y(dex); 
        [p,S] = polyfit(x,y,1); plot(x,p(1)*x+p(2),'m','LineWidth',3);
        x=optode_temperature;               y=optode_oxygen; 
        dex=find(y<sat100 & depth>15); x=x(dex); y=y(dex); 
        [p,S] = polyfit(x,y,1); plot(x,p(1)*x+p(2),'g','LineWidth',3);
        x=bottom.sea_water_temperature(ind);y=bottom.ctd_tc_oxygen(ind); 
        dex=find(y<sat100); x=x(dex); y=y(dex); 
        [p,S] = polyfit(x,y,1); plot(x,p(1)*x+p(2),'r','LineWidth',3);
        x=riser.sea_water_temperature(rnd); y=riser.ctd_tc_oxygen(rnd); 
        dex=find(y<sat100); x=x(dex); y=y(dex); 
        [p,S] = polyfit(x,y,1); plot(x,p(1)*x+p(2),'b','LineWidth',3);
        legend({'CSPP near NSIF','CSPP near MFN','NSIF','MFN','fit CSPP near NSIF','fit CSPP near MFN','fit NSIF','fit MFN'},'Location','northwest');
        ylim([0 sat100]); xlim([6.8 12]);
        set(gca,'fontsize',fsize); 
        title([int2str(max(profile)),' Profiles, Apr-Jun 2024'])
        
        set(gcf,'units','inches','position',[1 1 10 7]);

        % https://www.mathworks.com/help/stats/linearmodel.html
        x1=optode_temperature; x2=sea_water_practical_salinity; y=optode_oxygen; 
        dex=find(y<300); x1=x1(dex); x2=x2(dex); y=y(dex); 
        % mdl = fitlm([x1,x2],y)
        % beta = mvregress([x1,x2],y)
        [b,bint,r,rint,stats] = regress(y,[x2*0+1,x1,x2]); disp(b)
        % All methods give the same slopes and no intercept for large numbers of measurements. 
                
        close all;
        figure(103); 
        scatter3(x1,x2,y,'filled'); hold on
        x1fit = min(x1):.1:max(x1);
        x2fit = min(x2):.1:max(x2);
        [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
        YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT; % + b(4)*X1FIT.*X2FIT;
        mesh(X1FIT,X2FIT,YFIT);

        j=1; diffstats=NaN*zeros(max(profile),4);
        for j=1:max(profile)
            j=double(j);
            figure(j);
            dex=find(profile == j);
            if ~isempty(dex)
                lagtT=lagTime(convertTo(Time(dex),'posixtime'),tauT);
                lagfakeT=interp1(lagtT,optode_temperature(dex),convertTo(Time(dex),'posixtime'));
                lagT(dex)=lagfakeT;

                subplot(231);
                plot(sea_water_temperature(dex),depth(dex),'k',optode_temperature(dex),depth(dex),'r',lagfakeT,depth(dex),'g','LineWidth',2); grid on;
                xlabel('temperature');
                title('unlag method improves data from optode thermistor');
                set(gca,'ydir','rev'); ylim([0 25]);
                subplot(232);
                plot(sea_water_temperature(dex)-optode_temperature(dex),depth(dex),'r',sea_water_temperature(dex)-lagfakeT,depth(dex),'g'); grid on;
                xlabel('temperature difference'); legend({'T from CTD-optode','T from CTD-unlagged optode'})
                set(gca,'ydir','rev'); ylim([0 25]);

                % get MFN and NSIF values
                ind = find( (mean(Time(dex))-1/24/4) < bottom.Time & bottom.Time < (mean(Time(dex))+1/24/4) );
                rnd = find( (mean(Time(dex))-1/24/4) <  riser.Time &  riser.Time < (mean(Time(dex))+1/24/4) );

                subplot(234)
                plot(optode_oxygen(dex),depth(dex),'LineWidth',2); grid on; hold on;
                % ADD shifted oxygen
                % ADD MFN and NSIF
                plot(bottom.ctd_tc_oxygen(ind),24+ind*0,'xy','MarkerSize',14,'LineWidth',3);
                plot( riser.ctd_tc_oxygen(rnd), 7+rnd*0,'xb','MarkerSize',14,'LineWidth',3);
                set(gca,'ydir','reverse');

                subplot(235)
                scatter(optode_temperature(dex),optode_oxygen(dex),12,depth(dex),'filled');
                grid on; hold on;
                plot(bottom.sea_water_temperature(ind),bottom.ctd_tc_oxygen(ind),'xb','LineWidth',3);
                plot( riser.sea_water_temperature(rnd), riser.ctd_tc_oxygen(rnd),'xy','LineWidth',3);
                xlabel('optode temperature'); ylabel('oxygen (uM)'); colorbar;
                if ~isempty(ind); title(string(bottom.Time(ind(1)))); end

                % oxygen lag
                lagtO=lagTime(convertTo(Time(dex),'posixtime'),tauO);
                lagfakeO=interp1(lagtO,optode_oxygen(dex),convertTo(Time(dex),'posixtime'));
                lagO(dex)=lagfakeO;

                subplot(233);
                plot(optode_oxygen(dex),depth(dex),'k',lagfakeO,depth(dex),'g','linewidth',2); hold on;
                if ~isempty(ind); text(bottom.ctd_tc_oxygen(ind(1)),24,'MFN'); end
                if ~isempty(rnd); text( riser.ctd_tc_oxygen(rnd(1)), 7,'NSIF');end
                xlim([0 (max([max(optode_oxygen(dex)) max(bottom.ctd_tc_oxygen(ind)) max(riser.ctd_tc_oxygen(rnd))])+10)]); 
                ylim([0 25]);
                set(gca,'ydir','reverse');

                % compute stats. Linear fit
                x=sea_water_temperature(dex);
                y=optode_temperature(dex);
                [p,S] = polyfit(x,y,1);
                diffTstats(j,1)=S.rsquared;
                x=sea_water_temperature(dex);
                y=lagfakeT;
                dex=~isnan(x.*y);
                x=x(dex==1);
                y=y(dex==1);
                [p,S] = polyfit(x,y,1);
                diffTstats(j,2)=S.rsquared;
            end
        end
        %% show that tempeature lagging works over all profiles.
        figure
        lagT(lagT<1)=NaN;
        d=0:.25:24;
        for j=1:length(d)
            dex=find(d(j)<depth & depth<(d(j)+.25));
            mean_sea_water_temperature(j)=mean(sea_water_temperature(dex),'omitnan');
            mean_optode_temperature(j)=mean(optode_temperature(dex),'omitnan');
            mean_lagT(j)=mean(lagT(dex),'omitnan');
            mean_depth(j)=mean(depth(dex),'omitnan');
        end; mean_lagT(4)=NaN;
        subplot(121);
        %plot(sea_water_temperature,depth,'k.',optode_temperature,depth,'r.',lagT,depth,'g.'); hold on;
        plot(mean_sea_water_temperature,mean_depth,'k',mean_optode_temperature,mean_depth,'r',mean_lagT,mean_depth,'g','LineWidth',3); grid on;        
        set(gca,'fontsize',fsize);
        xlabel('temperature (Celsius)'); ylabel('depth (m)')
        title(['Average over ',int2str(max(profile)),' Profiles, Apr-Jun 2024'])
        set(gca,'ydir','rev');
        ylim([0 22.5]);
        legend({'CTD','optode','unlagged optode'},'Location','east');
        subplot(122);
        %plot(sea_water_temperature-optode_temperature,depth,'r.',sea_water_temperature-lagT,depth,'g.'); hold on;
        plot(mean_sea_water_temperature-mean_optode_temperature,mean_depth,'r',mean_sea_water_temperature-mean_lagT,mean_depth,'g','LineWidth',3); grid on;
        set(gca,'fontsize',16);
        xlabel('temperature difference'); legend({'T from CTD-optode T','T from CTD-unlagged optode T'},'Location','southeast');
        set(gca,'ydir','rev');
        ylim([0 22.5]); ylabel('depth (m)');
        title('Makes slow optode T more like fast CTD T');
        set(gcf,'units','inches','position',[1 1 12 7]);

        %% show r^2 improvement for temperature
        figure
        plot(diffTstats(:,1),diffTstats(:,2),'k.',[0.8 1],[0.8 1]); axis equal; grid on;
        xlabel('r^2 CTD T to optode T');
        ylabel('r^2 CTD T to unlagged optode T');
        title(num2str([median(diffTstats(:,1)) median(diffTstats(:,2))]))

        %% show that NSIF DO is better predicted for the first 76 profiles
        figure
        lagO(lagO<1)=NaN;
        d=0:.25:25;
        for j=1:length(d)
            dex=find(d(j)<depth & depth<(d(j)+.25) & profile<77);
            mean_optode_oxygen(j)=mean(optode_oxygen(dex),'omitnan');
            mean_lagO(j)=mean(lagO(dex),'omitnan');
            mean_depth(j)=mean(depth(dex),'omitnan');
        end 
        dex=find(profile<77);        % get MFN and NSIF values
        ind = find( min(Time(dex)) < bottom.Time & bottom.Time < max(Time(dex)) );
        rnd = find( min(Time(dex)) <  riser.Time &  riser.Time < max(Time(dex)) );
        plot(mean_optode_oxygen,mean_depth,'r',mean_lagO,mean_depth,'g','LineWidth',2); grid on; hold on;       
        text(median(bottom.ctd_tc_oxygen(ind)),24,'x MFN','FontSize',fsize); 
        text(median( riser.ctd_tc_oxygen(rnd)), 6,'x NSIF','FontSize',fsize);
        set(gca,'fontsize',fsize);
        xlabel('dissolved oxygen (uM)');
        title('Average over 76 Profiles, Apr-Jun 2024');
        set(gca,'ydir','rev');
        ylim([0 25.5]); ylabel('depth (m)'); xlim([100 360]);
        legend({'optode','unlagged optode'},'Location','southeast');
        set(gcf,'units','inches','position',[1 1 6 6]);
        
    end
end
toc

%% functions
function t1=lagTime(t,tau)
t1=t*0;
for i=1:length(t)
    denom=0;
    numer=0;
    for j=1:i
        numer=numer+t(j)*exp(-(t(i)-t(j))/tau);
        denom=denom+exp(-(t(i)-t(j))/tau);
    end
    t1(i)=numer/denom;
end
end

