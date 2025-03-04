clear; close all; clc;
tic
% https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/catalog.html
downloadTHREDDS=1;
addProfile=0;
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
            sensor = '10-PARADJ000';
            mooringSite = 'CE01ISSM';
            riserNode = 'RID16';
        elseif nsite == 2 % set up/update required inputs -- CE02SHSP
            site = 'CE02SHSP';
            sensor = '09-PARADJ000';
            mooringSite = 'CE02SHSM';
            riserNode= 'RID26';
        elseif nsite == 3    % set up/update required inputs -- CE06ISSP
            site = 'CE06ISSP';
            sensor = '10-PARADJ000';
            mooringSite = 'CE06ISSM';
            riserNode = 'RID16';
        else % set up/update required inputs -- CE07SHSP
            site = 'CE07SHSP';
            sensor = '09-PARADJ000';
            mooringSite = 'CE07SHSM';
            riserNode = 'RID26';
        end
        node = 'SP001';
        method = 'recovered_cspp';
        stream = 'parad_j_cspp_instrument_recovered';
        tag = '.*PARAD.*\.nc$';
        riserSensor =  '08-SPKIRB000';
        riserStream =  'spkir_abj_dcl_instrument_recovered';
        riserTag = '.*SPKIR.*\.nc$'; 
        mooringMethod = 'recovered_host';
        riser = load_gc_thredds_skip(mooringSite, riserNode, riserSensor, mooringMethod, riserStream, riserTag);
        riserannotations = get_annotations(mooringSite, riserNode, riserSensor);
        CSPPannotations =  get_annotations(site, node, sensor);
        disp('completed riser load');
        ce = load_gc_thredds_skip(site, node, sensor, method, stream, tag);
        disp('completed profiler load.');
        if nsite==1
            save('CE01ISSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        else
            save('CE07SHSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        end
        disp(int2str(nsite));
    end
end
disp(' Loaded THREDDS');

% ce.Properties.VariableNames' % view variable names

%% add profile variable to THREDDS
% no DOSTA data nsite= 1 deployment= 11
% no DOSTA data nsite= 1 deployment= 21
% no DOSTA data nsite= 2 deployment= 29
% no CTD data nsite= 2 deployment= 35
% no CTD data nsite= 3 deployment= 6
% no CTD data nsite= 4 deployment= 15
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
            save('CE01ISSPdosta.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSPdosta.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSPdosta.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        else
            save('CE07SHSPdosta.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        end
        close all;
    end
end
disp('  Added profile variable');

%% look at qc flags
% data site 1, deployment 24, bad calibration
% site 1, deployment 22. CE01ISSM-NSIF oxygen biofouling september 2023
% no data in site: 1 deployment: 21
% no data in site: 1 deployment: 11
% site 1, deployment 6. CE01ISSM-MFN oxygen biofouling aug sept 2016. 
% site 1, deployment 4. CE01ISSM-NSIF & MFN biofouling aug sept 2015.

% site 2, deployment 37. DO 10^6. 
% site 2, deployment 36. Too high. Slightly bad calibration.
% no data in site 2: deployment: 35
% site 2, deployment 34. bad calibration
% no data in site 2: deployment: 29
% site 2, deployment 8. NSIF biofouling.
% site 2, deployment 7. NSIF biofouling.
% site 2, deployment 5. NSIF biofouling.
% site 2, deployment 4. NSIF biofouling.

% site 3, deployment 12 July - September 2020 NSIF many zeros. Grass
% site 3, deployment 10 August 2019 NSIF zeroes out
% site 3, deployment 6 no data
% site 3, deployment 3 Aug Sept 2015 NSIF zeroes out. 
% site 3, deployment 2 Jul Aug 2015 NSIF zeroes out. 

% site 4, deployment 16. Bad calibration
% site 4, deployment 15 no data
% site 4, deployment 11 NSIF fields of grass July 2021
% site 4, deployment 2 NSIF biofoulinged August 2015

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
        for i=1:1:max(ce.deployment)
            figure(double(i));
            dex = find(ce.deployment==i);
            depth = ce.depth(dex);
            profiler_datetime = ce.profiler_datetime(dex);
            dissolved_oxygen = ce.dissolved_oxygen(dex);
            dosta_abcdjm_cspp_tc_oxygen=ce.dosta_abcdjm_cspp_tc_oxygen(dex);
            estimated_oxygen_concentration=ce.estimated_oxygen_concentration(dex);
            estimated_oxygen_saturation=ce.estimated_oxygen_saturation(dex);
            optode_temperature=ce.optode_temperature(dex);
            sea_water_temperature = ce.sea_water_temperature(dex);
            profile = ce.profile(dex);
            dissolved_oxygen_qartod_results = int16(ce.dissolved_oxygen_qartod_results(dex)); % flags least, but this is the one we'll use
            dosta_abcdjm_cspp_tc_oxygen_qartod_results = int16(ce.dosta_abcdjm_cspp_tc_oxygen_qartod_results(dex)); % flags most
            estimated_oxygen_concentration_qartod_results= int16(ce.estimated_oxygen_concentration_qartod_results(dex));
            % flags overlap, so most>middle>least
            % _qartod_results. look at values above 1.

            subplot(411);
            plot(profiler_datetime,dosta_abcdjm_cspp_tc_oxygen-estimated_oxygen_concentration,'.')
            hold on; plot(profiler_datetime,depth*0,'k',profiler_datetime,depth*0+.2,'k:',profiler_datetime,depth*0-.2,'k:'); grid on;
            title('oxygen (from instrument - calculated): should be near zero');
            % different by >.1 only in the last deployment. Why the 2 states? Bad calibration.
            % dissolved_oxygen vs. either not similar. This one is CTD corrected. Use this one.

            subplot(412);
            plot(profiler_datetime,sea_water_temperature-optode_temperature,'.');
            hold on; plot(profiler_datetime,depth*0,'k',profiler_datetime,depth*0+4,'r');  grid on;
            title('Temperature (CTD-optode): should average zero');
            % sea_water and optode temperatues are different.
            % Note if there is a problem with one of them.

            subplot(413);
            if ~isempty(dex)
                scatter(profiler_datetime,depth,10,dissolved_oxygen,'filled')
                ylabel('oxygen'); hold on; set(gca,'ydir','rev'); 
                a=colorbar; set(a,'limits',[0 400]);
                dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                if ~isempty(dex)
                    scatter(riser.Time(dex),dex*0+7,10,riser.dissolved_oxygen(dex));
                end
                dex=find(dissolved_oxygen_qartod_results>1);
                if ~isempty(dex)
                    plot(profiler_datetime(dex),dissolved_oxygen_qartod_results(dex),'r.','MarkerSize',12);
                end
                title('CTD corrected DO vs. moorings with red QARTOD');
            else
                disp(['no data in site: ',int2str(nsite),' deployment: ',int2str(i)]);
            end

            subplot(414);
            dex = find(ce.deployment==i);
            if ~isempty(dex)
                dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                plot(riser.Time(dex),riser.dissolved_oxygen(dex),'b-'); hold on;
                mintime=datetime(year(min(profiler_datetime)),month(min(profiler_datetime)),day(min(profiler_datetime)));
                maxtime=datetime(year(max(profiler_datetime)),month(max(profiler_datetime)),day(max(profiler_datetime)));
                dt=.25; % days
                tmptime=mintime:dt:(dt+maxtime);  tmptime.TimeZone = 'UTC';
                tmpoxygen=(1:length(tmptime))*NaN;
                for j=1:length(tmptime)
                    ind=find(tmptime(j)<=profiler_datetime & profiler_datetime<(tmptime(j)+dt)); % & depth < 5 & depth < 10);
                    if ~isempty(ind)
                        tmpoxygen(j)=mean(dissolved_oxygen(ind));
                    end
                end
                plot(tmptime,tmpoxygen,'go-'); grid on;
                if ~isempty(dex)
                    legend({'NSIF','CSPP'});
                end
                title('DO time series NSIF depth');
                ylim([0 400]);
            end

            set(gcf,'units','inches','position',[1 1 16 10]);
        end

        % if nsite==1
        %     save('CE01ISSPdosta.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % elseif nsite ==2
        %     save('CE02SHSPdosta.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % elseif nsite==3
        %     save('CE06ISSPdosta.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % else
        %     save('CE07SHSPdosta.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % end
    end
end

toc