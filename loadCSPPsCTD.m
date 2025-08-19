clear; close all; clc;
tic
% https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/catalog.html
downloadTHREDDS=1;
addProfile=1;
inspectQC=1;
folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder);
nsitedepths=[25,80,29,87];
sites=1:4;

%% gather THREDDS CSPP data
if downloadTHREDDS
    for nsite=sites
        if nsite==1        % set up required inputs -- CE01ISSP
            site = 'CE01ISSP';
            sensor = '09-CTDPFJ000';
            mooringSite = 'CE01ISSM';
            buoyNode='SBD17';
            riserNode = 'RID16';
            bottomNode = 'MFD37';
            buoySensor = '06-CTDBPC000';
            riserSensor = '03-CTDBPC000';
            bottomSensor = '03-CTDBPC000';
            buoyStream = 'ctdbp_cdef_instrument_recovered';
            riserStream = buoyStream; bottomStream = buoyStream;
            buoyTag = '.*CTDBP.*\.nc$';
        elseif nsite == 2 % set up/update required inputs -- CE02SHSP
            site = 'CE02SHSP';
            sensor = '08-CTDPFJ000';
            mooringSite = 'CE02SHSM';
            buoyNode = 'SBD11';
            riserNode= 'RID27';
            buoySensor = '06-METBKA001';
            riserSensor = '03-CTDBPC000';
            buoyStream = 'metbk_ct_instrument';
            riserStream = 'ctdbp_cdef_instrument_recovered';
            buoyTag = '.*METBK.*\.nc$';
            % CE02SHBP-LJ01D-06-CTDBPN106-streamed-ctdbp_no_sample. TOO DENSE. GET FROM BRANDY.
        elseif nsite == 3    % set up/update required inputs -- CE06ISSP
            site = 'CE06ISSP';
            sensor = '09-CTDPFJ000';
            mooringSite = 'CE06ISSM';
            buoyNode = 'SBD17';
            riserNode = 'RID16';
            bottomNode = 'MFD37';
            buoySensor = '06-CTDBPC000';
            riserSensor = '03-CTDBPC000';
            bottomSensor = '03-CTDBPC000';
            buoyStream = 'ctdbp_cdef_instrument_recovered';
            riserStream = buoyStream; bottomStream = buoyStream;
            buoyTag = '.*CTDBP.*\.nc$';
        else % set up/update required inputs -- CE07SHSP
            site = 'CE07SHSP';
            sensor = '08-CTDPFJ000';
            mooringSite = 'CE07SHSM';
            buoyNode = 'SBD11';
            riserNode = 'RID27';
            bottomNode = 'MFD37';
            buoySensor = '06-METBKA001';
            riserSensor = '03-CTDBPC000';
            bottomSensor = '03-CTDBPC000';
            buoyStream = 'metbk_ct_instrument';
            riserStream = 'ctdbp_cdef_instrument_recovered';
            bottomStream = riserStream;
            buoyTag = '.*METBK.*\.nc$';
        end
        node = 'SP001';
        method = 'recovered_cspp';
        stream = 'ctdpf_j_cspp_instrument_recovered';
        tag = '.*CTDPF.*\.nc$';
        riserTag = '.*CTDBP.*\.nc$'; bottomTag = riserTag;
        mooringMethod = 'recovered_inst';
        buoy = load_gc_thredds(mooringSite, buoyNode, buoySensor, mooringMethod, buoyStream, buoyTag);
        riser = load_gc_thredds(mooringSite, riserNode, riserSensor, mooringMethod, riserStream, riserTag);
        buoyannotations = get_annotations(mooringSite, buoyNode, buoySensor);
        riserannotations = get_annotations(mooringSite, riserNode, riserSensor);
        CSPPannotations =  get_annotations(site, node, sensor);

        if nsite==2
            load TSVel_NH10_2014_2024_V2;
            bottom=timetable(datetime(time,'ConvertFrom','datenum','TimeZone','UTC')',temp(41,:)',sal(41,:)');
            bottom.Properties.VariableNames{1}='sea_water_temperature';
            bottom.Properties.VariableNames{2}='sea_water_practical_salinity';
            bottomannotations = get_annotations('CE02SHBP', 'LJ01D','06-CTDBPN106');
        else
            bottom = load_gc_thredds(mooringSite, bottomNode, bottomSensor, mooringMethod, bottomStream, bottomTag);
            bottomannotations = get_annotations(mooringSite, bottomNode,bottomSensor);
        end

        ce = load_gc_thredds(site, node, sensor, method, stream, tag);
        % tag = 'deployment0017.*CTDPF.*\.nc$';

        if nsite==1
            save('CE01ISSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite ==2
            buoy.sea_surface_practical_salinity = gsw_SP_from_C ( buoy.sea_surface_conductivity*10, buoy.sea_surface_temperature, 1 );
            buoy.sea_surface_practical_salinity(buoy.sea_surface_practical_salinity<10)=NaN;
            save('CE02SHSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        else
            buoy.sea_surface_practical_salinity = gsw_SP_from_C ( buoy.sea_surface_conductivity*10, buoy.sea_surface_temperature, 1 );
            buoy.sea_surface_practical_salinity(buoy.sea_surface_practical_salinity<10)=NaN;
            save('CE07SHSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        end
        %filename = 'https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/CE02SHSP-SP001-08-CTDPFJ000-recovered_cspp-ctdpf_j_cspp_instrument_recovered.nc';
        %t = nc_reader(filename)
        %filesname = 'https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/deployment0001_CE02SHSP-SP001-08-CTDPFJ000-recovered_cspp-ctdpf_j_cspp_instrument_recovered_20150318T193656.769000-20150401T184856.125000.nc';
        %t = ncinfo(filename)
    end
end
disp(' Loaded THREDDS');

% ce.Properties.VariableNames' % view variable names

%% add profile variable to THREDDS

if 1==addProfile
    cd(folder);
    for nsite = sites
        if nsite==1
            load CE01ISSP.mat;
        elseif nsite ==2
            load CE02SHSP.mat;
        elseif nsite==3
            load CE06ISSP.mat;
        else
            load CE07SHSP.mat;
        end
        ce.profiler_datetime = datetime(ce.profiler_timestamp,'ConvertFrom','posixtime');
        ce.profiler_datetime.TimeZone = 'UTC';
        ce.profile = ce.deployment*NaN;
        % find profile number
        for i=1:max(ce.deployment)
            dex=find(i==ce.deployment);
            % only if length dex > 10. There are some data. None in ce06
            % D00006. 2024-08-23 check back for fix in December
            if length(dex)>10
                Time=ce.Time(dex);
                internal_timestamp=ce.internal_timestamp(dex);
                profiler_timestamp=ce.profiler_timestamp(dex);
                % figure('name',int2str(i));
                % plot(dex(2:end),diff(internal_timestamp),'bo',dex(2:end),diff(profiler_timestamp),'gx','MarkerSize',24); hold on;
                % plot(dex(2:end),diff(convertTo(Time,'posixtime')),'r.');
                % Time, internal_timestamp, and profiler_timestamp are the same data
                figure;
                tmp=diff(convertTo(Time,'posixtime'));
                scatter(dex(2:end),tmp,12,ce.sea_water_pressure(dex(2:end)),'filled');
                hold on;grid on;
                ind=find(tmp>300);
                ce.profile(dex(1):dex(ind(1)))=1;
                for j=2:length(ind)
                    plot([1 1]*dex(ind(j)),[0 300],'k');
                    ce.profile(dex((1+ind(j-1)):ind(j)))=j;
                end
                ylim([0 1/8]);
                ce.profile(dex((ind(j)+1):end))=length(ind)+1;
                figure
                plot(ce.Time(dex),ce.profile(dex),'.');
                % looks like diff 5 minutes is a good cut off.
                close all;
            else
                disp(['no data nsite= ',int2str(nsite),' deployment= ',int2str(i)]);
                % no data nsite= 2 deployment= 35
                % no data nsite= 3 deployment= 6
                % no data nsite= 4 deployment= 15
            end
        end
        if nsite==1
            save('CE01ISSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        else
            save('CE07SHSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        end
    end
end
disp('  Added profile variable');

%% look at qc flags
% dex=find(ce.sea_water_temperature_qartod_results>1 | ce.sea_water_practical_salinity_qartod_results>1);
% figure(777);
% scatter(ce.profiler_datetime(dex),ce.depth(dex),10,ce.sea_water_practical_salinity(dex),'filled');
% colorbar;
% title('Chris has fixes that CI is reviewing. Come back to this later.');
% 2024-08-23 flags from ce07 not showing up yet. Updates are in the batch CI is testing.
close all;
% nsite 1 bad surface salinity 4360
% nsite 2 bad depth 1
% nsite 2 bad salinity 2387
% nsite 2 bad surface salinity 2564
% nsite 2 no data deployment 35
% nsite 3 bad salinity 39727
% nsite 3 bad surface salinity 42189
% nsite 3 no data deployment 6
% nsite 4 bad salinity 72
% nsite 4 no data deployment 15
if inspectQC
    cd(folder);
    for nsite = sites
        if nsite==1
            load CE01ISSP.mat;
        elseif nsite ==2
            load CE02SHSP.mat;
        elseif nsite==3
            load CE06ISSP.mat;
        else
            load CE07SHSP.mat;
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
            if isempty(dex)
                disp(['nsite ',int2str(nsite),' no data deployment ',int2str(i)]);
            else
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
                    if ~(nsite == 3 && i == 6) %% this one bad deployment will be reloaded
                        subplot(221);
                        scatter(profiler_datetime,depth,10,sea_water_temperature,'filled');
                        ylabel('temperature'); hold on; set(gca,'ydir','rev'); colorbar;
                        dex=find(min(profiler_datetime)<buoy.Time & buoy.Time<max(profiler_datetime));
                        scatter(buoy.Time(dex),dex*0+1,10,buoy.sea_water_temperature(dex));
                        dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                        scatter(riser.Time(dex),dex*0+7,10,riser.sea_water_temperature(dex));
                        dex=find(min(profiler_datetime)<bottom.Time & bottom.Time<max(profiler_datetime));
                        scatter(bottom.Time(dex),dex*0+nsitedepths(nsite)-1,10,bottom.sea_water_temperature(dex));
                        subplot(222);
                        ylabel('salinity'); hold on; set(gca,'ydir','rev'); colorbar;
                        scatter(profiler_datetime,depth,10,sea_water_practical_salinity,'filled');
                        dex=find(min(profiler_datetime)<buoy.Time & buoy.Time<max(profiler_datetime));
                        scatter(buoy.Time(dex),dex*0+1,10,buoy.sea_water_practical_salinity(dex));
                        dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                        scatter(riser.Time(dex),dex*0+7,10,riser.sea_water_practical_salinity(dex));
                        dex=find(min(profiler_datetime)<bottom.Time & bottom.Time<max(profiler_datetime));
                        scatter(bottom.Time(dex),dex*0+nsitedepths(nsite)-1,10,bottom.sea_water_practical_salinity(dex));

                        % 1:1 plots
                        for p=1:max(profile)
                            dex=find(1<depth & depth<2 & p==profile);
                            csppT.buoy(nsite,i,p)=mean(sea_water_temperature(dex));
                            ind=find(p==profile);
                            dex=find((min(profiler_datetime(ind))-15/60/24)<buoy.Time & buoy.Time<(max(profiler_datetime(ind))+15/60/24));
                            moorT.buoy(nsite,i,p)=mean(buoy.sea_water_temperature(dex));

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
                        plot(squeeze(csppT.buoy(nsite,i,:)),squeeze(moorT.buoy(nsite,i,:)),'.',squeeze(csppT.riser(nsite,i,:)),squeeze(moorT.riser(nsite,i,:)),'.',squeeze(csppT.bottom(nsite,i,:)),squeeze(moorT.bottom(nsite,i,:)),'.');
                        grid on; axis equal;
                    end
                else % 2 or 4
                    subplot(211)
                    scatter(profiler_datetime,depth,10,sea_water_temperature,'filled');
                    ylabel('temperature'); hold on; set(gca,'ydir','rev'); colorbar;
                    dex=find(min(profiler_datetime)<buoy.Time & buoy.Time<max(profiler_datetime));
                    scatter(buoy.Time(dex),dex*0+1,10,buoy.sea_surface_temperature(dex));
                    dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                    scatter(riser.Time(dex),dex*0+7,10,riser.sea_water_temperature(dex));
                    dex=find(min(profiler_datetime)<bottom.Time & bottom.Time<max(profiler_datetime));
                    scatter(bottom.Time(dex),dex*0+nsitedepths(nsite)-1,10,bottom.sea_water_temperature(dex));
                    subplot(212);
                    scatter(profiler_datetime,depth,10,sea_water_practical_salinity,'filled')
                    ylabel('salinity'); hold on; set(gca,'ydir','rev'); colorbar;
                    dex=find(min(profiler_datetime)<buoy.Time & buoy.Time<max(profiler_datetime));
                    scatter(buoy.Time(dex),dex*0+1,10,buoy.sea_surface_practical_salinity(dex));
                    dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                    scatter(riser.Time(dex),dex*0+7,10,riser.sea_water_practical_salinity(dex));
                    dex=find(min(profiler_datetime)<bottom.Time & bottom.Time<max(profiler_datetime));
                    scatter(bottom.Time(dex),dex*0+nsitedepths(nsite)-1,10,bottom.sea_water_practical_salinity(dex));
                end
                close; % takes too much memory to plot them all.
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
            save('CE01ISSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        else
            save('CE07SHSP.mat','-v7.3','ce','buoy','riser','bottom','buoyannotations','riserannotations','bottomannotations','CSPPannotations');
        end
    end
end

toc
% filter: Temperature 0.085 sec, Conductivity 0.085 sec, Pressure 0.25 sec.
% align
%   temperature advance = 0.0625 seconds (1/16). C is instant measurements with 1.04 cm 0.04 sec spatial delay, P is instant.
%   (TAadvance not part of the protocol)
% thermal mass
%   celltm Alpha = 0.03
%   celltm Tau (1/beta)= 7.0
% skip loop edit because this profiler is different
% calculate salinity, density, etc
% bin average

