clear; close all; clc;
tic

folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder);
% ce.Properties.VariableNames' % view variable names

load CE01ISSPdosta.mat;
ctd=load('CE01ISSP');


bottomData=timetable(datetime(ce.profiler_timestamp(1),'ConvertFrom','datenum','TimeZone','UTC')',0,0,0,0,0);
bottomData.Properties.VariableNames{1}='profiler_timestamp';
bottomData.Properties.VariableNames{2}='sea_water_practical_salinity';
bottomData.Properties.VariableNames{3}='sea_water_temperature';
bottomData.Properties.VariableNames{4}='dissolved_oxygen';
bottomData.Properties.VariableNames{5}='maxdepth';

%bottomData.profiler_timestamp=zeros(1,numprofiles)*NaN;
%bottomData.sea_water_practical_salinity=bottomData.profiler_timestamp;
%bottomData.sea_water_temperature=bottomData.profiler_timestamp;
%bottomData.dissolved_oxygen=bottomData.profiler_timestamp;
%bottomData.maxdepth =bottomData.profiler_timestamp;
cnt=1;

for i=1:max(ctd.ce.deployment)
    dex=find(i==ctd.ce.deployment);
    if ~isempty(dex)
        ctd_Time=ctd.ce.Time(dex);
        ctd_profiler_timestamp=ctd.ce.profiler_timestamp(dex);
        ctd_profile = ctd.ce.profile(dex);
        ctd_sea_water_practical_salinity = ctd.ce.sea_water_practical_salinity(dex);
        ctd_sea_water_temperature = ctd.ce.sea_water_temperature(dex);
        ctd_depth = ctd.ce.depth(dex);
        % Time, internal_timestamp, and profiler_timestamp are the same data
        dex = find(i==ce.deployment);
        if ~isempty(dex)
            for j=1:max(ctd_profile)
                ind = find(j == ctd_profile);
                if ~isempty(ind)
                    index = find(ctd_depth(ind)<(max(ctd_depth(ind))-2));
                    bottomData.profiler_timestamp(cnt)= mean(ctd_profiler_timestamp(ind(index)));
                    bottomData.sea_water_practical_salinity(cnt)= mean(ctd_sea_water_practical_salinity(ind(index)));
                    bottomData.sea_water_temperature(cnt)= mean(ctd_sea_water_temperature(ind(index)));
                    bottomData.maxdepth(cnt)= max(ctd_depth(ind));
                    bottomData.Time(cnt)= min(ctd_Time(ind));
                    bottomData.profiler_timestamp(cnt)= min(ctd_profiler_timestamp(ind));

                    ind = find(j == ce.profile(dex));
                    if ~isempty(ind)
                        index = find(ce.depth(dex(ind))<(bottomData.maxdepth(cnt)-2));
                        bottomData.dissolved_oxygen(cnt)=mean(ce.dissolved_oxygen(dex(ind(index))));
                    end
                    cnt=cnt+1;
                end
            end
        else
            disp(['no DOSTA data deployment = ',int2str(i)]);
        end
    else
        disp(['no CTD data deployment = ',int2str(i)]);
    end
    disp(int2str(i));
end
disp(int2str(cnt));

figure
subplot(221); plot(bottomData.Time,bottomData.sea_water_practical_salinity,'k.');
ylabel('salinity');
subplot(222); plot(bottomData.Time,bottomData.sea_water_temperature,'k.');
ylabel('temperature')
subplot(223); plot(bottomData.Time,bottomData.dissolved_oxygen,'k.');
ylabel('dissolved oxygen');
subplot(224); plot(bottomData.Time,bottomData.maxdepth,'k.');
ylabel('max profile depth');

save('CE01ISSPbottomData.mat','-v7.3','bottomData');

toc