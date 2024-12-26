clear; close all; clc;
tic
% https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/catalog.html
downloadTHREDDS=0;
addProfile=1;
inspectQC=0;
folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder);
nsitedepths=[25,80,29,87];
sites=1:4;

%% gather THREDDS CSPP data
if downloadTHREDDS
    for nsite=sites
        if nsite==1        % set up required inputs -- CE01ISSP
            site = 'CE01ISSP';
            sensor = '02-DOSTAJ000';
            mooringSite = 'CE01ISSM';
            riserNode = 'RID16';
            bottomNode = 'MFD37';
            riserSensor =  '03-DOSTAD000';
            bottomSensor = '03-DOSTAD000';
            riserStream =  'dosta_abcdjm_ctdbp_instrument_recovered';
            bottomStream = riserStream;
        elseif nsite == 2 % set up/update required inputs -- CE02SHSP
            site = 'CE02SHSP';
            sensor = '01-DOSTAJ000';
            mooringSite = 'CE02SHSM';
            riserNode= 'RID27';
            riserSensor = '04-DOSTAD000';
            riserStream = 'dosta_abcdjm_dcl_instrument_recovered';
            % CE02SHBP-LJ01D-06-CTDBPN106-streamed-ctdbp_no_sample. TOO DENSE. GET FROM BRANDY.
        elseif nsite == 3    % set up/update required inputs -- CE06ISSP
            site = 'CE06ISSP';
            sensor = '02-DOSTAJ000';
            mooringSite = 'CE06ISSM';
            riserNode = 'RID16';
            bottomNode = 'MFD37';
            riserSensor =  '03-DOSTAD000';
            bottomSensor = '03-DOSTAD000';
            riserStream = 'dosta_abcdjm_ctdbp_instrument_recovered';
            bottomStream = riserStream;
        else % set up/update required inputs -- CE07SHSP
            site = 'CE07SHSP';
            sensor = '01-DOSTAJ000';
            mooringSite = 'CE07SHSM';
            riserNode = 'RID27';
            bottomNode = 'MFD37';
            riserSensor = '04-DOSTAD000';
            bottomSensor = '03-DOSTAD000';
            riserStream = 'dosta_abcdjm_dcl_instrument_recovered';
            bottomStream ='dosta_abcdjm_ctdbp_instrument_recovered';
        end
        node = 'SP001';
        method = 'recovered_cspp';
        stream = 'dosta_abcdjm_cspp_instrument_recovered';
        tag = '.*DOSTA.*\.nc$';
        riserTag = '.*DOSTA.*\.nc$'; bottomTag = riserTag;
        if nsite == 1
            mooringMethod = 'recovered_inst';
        elseif nsite==2
            mooringMethod = 'recovered_host';
        elseif nsite == 3
            mooringMethod = 'recovered_inst';
        else
            mooringMethod = 'recovered_host';
        end
        riser = load_gc_thredds_skip(mooringSite, riserNode, riserSensor, mooringMethod, riserStream, riserTag);
        riserannotations = get_annotations(mooringSite, riserNode, riserSensor);
        CSPPannotations =  get_annotations(site, node, sensor);

        if nsite == 4
            mooringMethod = 'recovered_inst';
        end
        if nsite==2
            % MAKE A BOTTOM DOSTA FILE
            % load TSVel_NH10_2014_2024_V2;
            % bottom=timetable(datetime(time,'ConvertFrom','datenum','TimeZone','UTC')',temp(41,:)',sal(41,:)');
            % bottom.Properties.VariableNames{1}='sea_water_temperature';
            % bottom.Properties.VariableNames{2}='sea_water_practical_salinity';
            bottom = [];
            bottomannotations = get_annotations('CE02SHBP', 'LJ01D','06-DOSTAD106');
        else
            bottom = load_gc_thredds_skip(mooringSite, bottomNode, bottomSensor, mooringMethod, bottomStream, bottomTag);
            bottomannotations = get_annotations(mooringSite, bottomNode,bottomSensor);
        end
        ce = load_gc_thredds_skip(site, node, sensor, method, stream, tag);

        if nsite==1
            save('CE01ISSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        else
            save('CE07SHSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        end
    end
end
disp(' Loaded THREDDS');

% ce.Properties.VariableNames' % view variable names

%% add profile variable to THREDDS

if 1==addProfile
    cd(folder);
    for nsite = sites
        if nsite==1
            load CE01ISSPdosta.mat;
            ctd=load('CE01ISSP');
        elseif nsite ==2
            load CE02SHSPdosta.mat;
            ctd=load('CE02SHSP');
        elseif nsite==3
            load CE06ISSPdosta.mat;
            ctd=load('CE06ISSP');
        else
            load CE07SHSPdosta.mat;
            ctd=load('CE07SHSP');
        end
        ce.profiler_datetime = datetime(ce.profiler_timestamp,'ConvertFrom','posixtime');
        ce.profiler_datetime.TimeZone = 'UTC';
        ce.profile = ce.deployment*NaN;
        % find profile number
        for i=1:max(ctd.ce.deployment)
            dex=find(i==ctd.ce.deployment);
            if ~isempty(dex)
                ctd_Time=ctd.ce.Time(dex);
                ctd_profiler_timestamp=ctd.ce.profiler_timestamp(dex);
                ctd_profile = ctd.ce.profile(dex);
                % Time, internal_timestamp, and profiler_timestamp are the same data
                dex = find(i==ce.deployment);
                Time=ctd.ce.Time(dex);
                profiler_timestamp=ce.profiler_timestamp(dex);
                if ~isempty(dex)                    
                    for j=1:max(ctd_profile)
                        ind = find(j == ctd_profile);
                        index = find(min(ctd_profiler_timestamp(ind)) <= profiler_timestamp & profiler_timestamp <= max(ctd_profiler_timestamp(ind)));
                        ce.profile(dex(index))=j;
                    end
                    figure(double(i));
                    plot(ce.Time(dex),ce.profile(dex),'k.');
                    ylabel('profile number'); grid on; hold on;
                    title(int2str(i));
                    set(gcf,'units','inches','Position',[(1+double(i)/25) 1 12 8]);
                    plot(ctd_Time,ctd_profile,'ro');
                    legend({'DOSTA','CTD'},'Location','east');
                else
                    disp(['no DOSTA data nsite= ',int2str(nsite),' deployment= ',int2str(i)]);
                end
            else
                disp(['no CTD data nsite= ',int2str(nsite),' deployment= ',int2str(i)]);
            end
        end
        if nsite==1
            save('CE01ISSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        else
            save('CE07SHSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        end
    end
end
disp('  Added profile variable');

%% look at qc flags

if inspectQC
    cd(folder);
    for nsite = sites
        if nsite==1
            load CE01ISSPdosta.mat;
            ctd=load('CE01ISSP');
        elseif nsite ==2
            load CE02SHSPdosta.mat;
            ctd=load('CE02SHSP');
        elseif nsite==3
            load CE06ISSPdosta.mat;
            ctd=load('CE06ISSP');
        else
            load CE07SHSPdosta.mat;
            ctd=load('CE07SHSP');
        end
        dex=find(ce.depth>90);
        if ~isnan(dex)
            ce.depth(dex)=NaN;
            % nsite 2, deployment 4 has a bad depth near the profile bottom. Should be 69.9;
            disp(['nsite ',int2str(nsite),' bad depth ',int2str(length(dex))]);
        end
        if nsite==3
            dex=find(ce.sea_water_practical_salinity<19 | 35<ce.sea_water_practical_salinity);
        else
            dex=find(ce.sea_water_practical_salinity<21 | 35<ce.sea_water_practical_salinity);
        end
        if ~isnan(dex)
            ce.sea_water_practical_salinity(dex)=NaN;
            disp(['nsite ',int2str(nsite),' bad salinity ',int2str(length(dex))]);
        end
        if nsite<4
            dex=find(1<ce.sea_water_practical_salinity_qartod_results & ce.depth<1);
            if ~isnan(dex)
                ce.sea_water_practical_salinity(dex)=NaN;
                disp(['nsite ',int2str(nsite),' bad surface salinity ',int2str(length(dex))]);
            end
        end
        for i=1:1:max(ce.deployment)
            figure
            % for j=1:max(ce.profile(dex)) % is this necessary? Just group by deployment.
            %dex = find(ce.deployment==i & ce.profile==j);
            dex = find(ce.deployment==i);
            depth = ce.depth(dex);
            profiler_datetime = ce.profiler_datetime(dex);
            sea_water_density = ce.sea_water_density(dex);
            sea_water_electrical_conductivity = ce.sea_water_electrical_conductivity(dex);
            % sea_water_pressure = ce.sea_water_pressure(dex); % redundant
            sea_water_practical_salinity = ce.sea_water_practical_salinity(dex);
            sea_water_temperature = ce.sea_water_temperature(dex);
            profile = ce.profile(dex);
            if nsite < 4
                sea_water_temperature_qartod_results = int16(ce.sea_water_temperature_qartod_results(dex));
                sea_water_practical_salinity_qartod_results = int16(ce.sea_water_practical_salinity_qartod_results(dex));
            end
            %% look at qc. see above.
            % _qc_executed deprecated
            % _qc_results deprecated
            % _qartod_executed. ignore
            % _qartod_results. look at values above 1.
            %  all within 2m of surface? Not 2024-08-23. Check back later.

            if nsite == 1 || nsite == 3
                if ~(nsite == 3 && i == 6)
                    subplot(221);
                    scatter(profiler_datetime,depth,10,sea_water_temperature,'filled');
                    ylabel('temperature'); hold on; set(gca,'ydir','rev'); colorbar;
                    dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                    scatter(riser.Time(dex),dex*0+7,10,riser.sea_water_temperature(dex));
                    dex=find(min(profiler_datetime)<bottom.Time & bottom.Time<max(profiler_datetime));
                    scatter(bottom.Time(dex),dex*0+nsitedepths(nsite)-1,10,bottom.sea_water_temperature(dex));
                    subplot(222);
                    ylabel('salinity'); hold on; set(gca,'ydir','rev'); colorbar;
                    scatter(profiler_datetime,depth,10,sea_water_practical_salinity,'filled');
                    dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                    scatter(riser.Time(dex),dex*0+7,10,riser.sea_water_practical_salinity(dex));
                    dex=find(min(profiler_datetime)<bottom.Time & bottom.Time<max(profiler_datetime));
                    scatter(bottom.Time(dex),dex*0+nsitedepths(nsite)-1,10,bottom.sea_water_practical_salinity(dex));

                    % 1:1 plots
                    for p=1:max(profile)
                        dex=find(6<depth & depth<8 & p==profile);
                        csppT.riser(nsite,i,p)=mean(sea_water_temperature(dex));
                        ind=find(p==profile);
                        dex=find((min(profiler_datetime(ind))-15/60/24)<riser.Time & riser.Time<(max(profiler_datetime(ind))+15/60/24));
                        moorT.riser(nsite,i,p)=mean(riser.sea_water_temperature(dex));

                        dex=find((nsitedepths(nsite)-10)<depth & p==profile);
                        csppT.bottom(nsite,i,p)=mean(sea_water_temperature(dex));
                        ind=find(p==profile);
                        dex=find((min(profiler_datetime(ind))-15/60/24)<bottom.Time & bottom.Time<(max(profiler_datetime(ind))+15/60/24));
                        moorT.bottom(nsite,i,p)=mean(bottom.sea_water_temperature(dex));
                    end
                    subplot(223)
                    plot(squeeze(csppT.riser(nsite,i,:)),squeeze(moorT.riser(nsite,i,:)),'.',squeeze(csppT.bottom(nsite,i,:)),squeeze(moorT.bottom(nsite,i,:)),'.');
                    grid on; axis equal;
                end
            else % 2 or 4
                subplot(211)
                scatter(profiler_datetime,depth,10,sea_water_temperature,'filled');
                ylabel('temperature'); hold on; set(gca,'ydir','rev'); colorbar;
                dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                scatter(riser.Time(dex),dex*0+7,10,riser.sea_water_temperature(dex));
                dex=find(min(profiler_datetime)<bottom.Time & bottom.Time<max(profiler_datetime));
                scatter(bottom.Time(dex),dex*0+nsitedepths(nsite)-1,10,bottom.sea_water_temperature(dex));
                subplot(212);
                scatter(profiler_datetime,depth,10,sea_water_practical_salinity,'filled')
                ylabel('salinity'); hold on; set(gca,'ydir','rev'); colorbar;
                dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                scatter(riser.Time(dex),dex*0+7,10,riser.sea_water_practical_salinity(dex));
                dex=find(min(profiler_datetime)<bottom.Time & bottom.Time<max(profiler_datetime));
                scatter(bottom.Time(dex),dex*0+nsitedepths(nsite)-1,10,bottom.sea_water_practical_salinity(dex));
            end
        end
        % subplot(223);
        % ylabel('temperature qartod executed'); hold on; set(gca,'ydir','rev'); colorbar;
        % scatter(profiler_datetime,depth,10,sea_water_temperature_qartod_results,'filled');
        % subplot(224);
        % ylabel('salinity qartod results'); hold on; set(gca,'ydir','rev'); colorbar;
        % scatter(profiler_datetime,depth,10,sea_water_practical_salinity_qartod_results,'filled');
        % Salinity range 21 - 35 is a good mask, not the qartod output.
        % Salinity goes down to 19 June 15-18, 2017. Seen in nsite 3
        % mask within 1m, if salinity qartod>1 and salinity <25
        % Ignore the temperature qartod. It flags too much. 2023-09-23 Chris has wider fail ranges in the queue
        % No need to mask temperature range. It never fails.
        % nsite 3, deployment 5 CTD must have been higher because every profile has a ~10 cm low salinity zone
        % nsite 3, deployment 9 and 14 have a salinity streak that should be masked.
        % nsite 2, deployment 13 has a bad streak. Ignore single bad streaks. Annotations are supposed to be issues that last longer.
        % nsite 1, deployment 18 has a few low salinity profiles. Leave. Probably real.
        % nsite 1, deployment 4 has many bad streaks. perhaps clog
        % At the top of the water colunn, temperatures are fine even when salinities are not fine. I guess, sucking in air

        %% confirmed depth vs. sea_water_pressure
        % plot(depth,sea_water_pressure-depth,'k.'); % --> zero
        % if j==1; ylabel('pressure'); xlabel('pressure-depth'); axis equal; grid on; hold on; end

        %% Note to users CTD bottom skipping in ce06: 1, 5, 10, 16.
        % ce07 all good. ce01 4-19 not great. ce02: 10
        % These are times when it didn't log. It's real.

        % h=plot_sticks(t,z_c-dz,salt_c,.25);
        set(gcf,'units','inches','position',[1 1 16 10]);
        if nsite==1
            save('CE01ISSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        else
            save('CE07SHSPdosta.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
        end
    end
end

toc