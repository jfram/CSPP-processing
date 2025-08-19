clear; close all; clc;
tic

% Lisa June 26, Aug 4

% test method

dz=.25;
depth=20; % d and thus t go to infinity, lag goes to tau. 
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

disp('stop for now');

% do for CE01 and CE02
% plot TS for mooring and profiler around a profile
% plot T-DO. Note problem. Use T from DOSTA vs. T from CTD to show issue.
% make depth of DO represent time lag weight

folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder);
% ce.Properties.VariableNames' % view variable names
load CE01ISSPdosta.mat;

for i=1:max(ce.deployment)
    % grab DO 
    % plot DO in time for the whole deployment for all to see if reasonable
    % grab hour before through after
    % 131: plot T:S
    % 132: plot T:DO
    % rework DO
    %  for each time, associate weighted DO_weighted time
    %  Note: move the weighted DO back to the weighted time
    % 133: replot T:DO

end

toc

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

