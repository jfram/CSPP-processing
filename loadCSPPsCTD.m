clear; close all; clc;
tic
% https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/catalog.html
loadTHREDDS=1;
loadProfile=0;
cd('C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc');

%% gather THREDDS CSPP data
if loadTHREDDS
    for nsite=1:4 %4
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
        if nsite==2
            bottom = [];
        else
            bottom = load_gc_thredds(mooringSite, bottomNode, bottomSensor, mooringMethod, bottomStream, bottomTag);
        end
        ce = load_gc_thredds(site, node, sensor, method, stream, tag);
        % tag = 'deployment0017.*CTDPF.*\.nc$';

        if nsite==1
            save('CE01ISSP.mat','-v7.3','ce','buoy','riser','bottom');
        elseif nsite ==2
            save('CE02SHSP.mat','-v7.3','ce','buoy','riser','bottom');
        elseif nsite==3
            save('CE06ISSP.mat','-v7.3','ce','buoy','riser','bottom');
        else
            save('CE07SHSP.mat','-v7.3','ce','buoy','riser','bottom');
        end
        %filename = 'https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/CE02SHSP-SP001-08-CTDPFJ000-recovered_cspp-ctdpf_j_cspp_instrument_recovered.nc';
        %t = nc_reader(filename)
        %filesname = 'https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/deployment0001_CE02SHSP-SP001-08-CTDPFJ000-recovered_cspp-ctdpf_j_cspp_instrument_recovered_20150318T193656.769000-20150401T184856.125000.nc';
        %t = ncinfo(filename)
    end
end
disp(' loaded THREDDS');

%% view
ce.Properties.VariableNames'

%% add profile variable to THREDDS
for nsite = 1:4
    cd('C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB');
    if nsite==1
        load CE01ISSP.mat;
    elseif nsite ==2
        load CE02SHSP.mat;
    elseif nsite==3
        load CE06ISSP.mat;
    else
        load CE07SHSP.mat;
    end
    if 0==loadProfile
        ce.profiler_datetime = datetime(ce.profiler_timestamp,'ConvertFrom','posixtime');
        ce.profile = ce.deployment*NaN;
        % find profile number
        for i=1:max(ce.deployment)
            dex=find(i==ce.deployment);
            % only if length dex > 10. There are some data. None in ce06 D00006
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
            end
        end
        if nsite==1
            save('CE01ISSP.mat','-v7.3','ce','buoy','riser','bottom');
        elseif nsite ==2
            save('CE02SHSP.mat','-v7.3','ce','buoy','riser','bottom');
        elseif nsite==3
            save('CE06ISSP.mat','-v7.3','ce','buoy','riser','bottom');
        else
            save('CE07SHSP.mat','-v7.3','ce','buoy','riser','bottom');
        end
    end
end


%% look at qc flags
% dex=find(ce.sea_water_temperature_qartod_results>1 | ce.sea_water_practical_salinity_qartod_results>1);
% figure(777);
% scatter(ce.profiler_datetime(dex),ce.depth(dex),10,ce.sea_water_practical_salinity(dex),'filled');
% colorbar;
% title('Chris has fixes that CI is reviewing. Come back to this later.');
% 2024-08-23 flags from ce07 not showing up yet. Updates are in the batch CI is testing.

for i=1:1:max(ce.deployment)
    figure
    dex=find(ce.deployment == i);
    for j=1:max(ce.profile(dex))
        dex = find(ce.deployment==i & ce.profile==j);
        depth = ce.depth(dex);
        profiler_timestamp = ce.profiler_timestamp(dex);
        profiler_datetime = ce.profiler_datetime(dex);

        sea_water_density = ce.sea_water_density(dex);

        sea_water_electrical_conductivity = ce.sea_water_electrical_conductivity(dex);
        % sea_water_pressure = ce.sea_water_pressure(dex); % redundant
        sea_water_practical_salinity = ce.sea_water_practical_salinity(dex);
        sea_water_temperature = ce.sea_water_temperature(dex);
        sea_water_temperature_qartod_results = int16(ce.sea_water_temperature_qartod_results(dex));
        sea_water_practical_salinity_qartod_results = int16(ce.sea_water_practical_salinity_qartod_results(dex));
        %% look at qc. see above.
        % _qc_executed deprecated
        % _qc_results deprecated
        % _qartod_executed. ignore
        % _qartod_results. look at values above 1. 
        %  No qartod results from CE07. Check back later.
        % all within 2m of surface? Not 2024-08-23. Check back later.

        %%
        subplot(221);
        if j==1; ylabel('temperature'); hold on; set(gca,'ydir','rev'); colorbar; end
        scatter(profiler_datetime,depth,10,sea_water_temperature,'filled');
        subplot(222);
        if j==1; ylabel('salinity'); hold on; set(gca,'ydir','rev'); colorbar; end
        scatter(profiler_datetime,depth,10,sea_water_practical_salinity,'filled')
        subplot(223);
        if j==1; ylabel('temperature qartod executed'); hold on; set(gca,'ydir','rev'); colorbar; end
        scatter(profiler_datetime,depth,10,sea_water_temperature_qartod_results,'filled');
        subplot(224);
        if j==1; ylabel('salinity qartod results'); hold on; set(gca,'ydir','rev'); colorbar; end
        scatter(profiler_datetime,depth,10,sea_water_practical_salinity_qartod_results,'filled');

        %% confirmed depth vs. sea_water_pressure
        % plot(depth,sea_water_pressure-depth,'k.'); % --> zero
        % if j==1; ylabel('pressure'); xlabel('pressure-depth'); axis equal; grid on; hold on; end

        %% Note to users CTD bottom skipping in ce06: 1, 5, 10, 16.
        % ce07 all good. ce01 4-19 not great. ce02: 10
        % bad depth in ce02 deployment 4. These are times when it didn't
        % log. It's real.

        %% save each deployment to a structure
        % profs{j}.each_parameter

    end
    % loop again to plot
    % h=plot_sticks(t,z_c-dz,salt_c,.25);
    set(gcf,'units','inches','position',[1 1 16 10]);
end

%% plot sticks for one deployment
% make a structure for storing profile

toc