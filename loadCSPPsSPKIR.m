clear; close all; clc;
tic
% 412, 443, 490, 510, 555, 619, 684 
% https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/catalog.html
downloadTHREDDS=0;
addProfile=0;
inspectQC=1;
folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder);
nsitedepths=[25,80,29,87];
sites=1; %1:4;

%% gather THREDDS CSPP data
if downloadTHREDDS
    for nsite=sites
        if nsite==1        % set up required inputs -- CE01ISSP
            site = 'CE01ISSP';
            sensor = '07-SPKIRJ000';
            mooringSite = 'CE01ISSM';
            riserNode = 'RID16';
        elseif nsite == 2 % set up/update required inputs -- CE02SHSP
            site = 'CE02SHSP';
            sensor = '06-SPKIRJ000';
            mooringSite = 'CE02SHSM';
            riserNode= 'RID26';
        elseif nsite == 3    % set up/update required inputs -- CE06ISSP
            site = 'CE06ISSP';
            sensor = '07-SPKIRJ000';
            mooringSite = 'CE06ISSM';
            riserNode = 'RID16';
        else % set up/update required inputs -- CE07SHSP
            site = 'CE07SHSP';
            sensor = '06-SPKIRJ000';
            mooringSite = 'CE07SHSM';
            riserNode = 'RID26';
        end
        node = 'SP001';
        method = 'recovered_cspp';
        stream = 'spkir_abj_cspp_instrument_recovered';
        tag = '.*SPKIR.*\.nc$';
        riserSensor =  '08-SPKIRB000';
        riserStream =  'spkir_abj_dcl_instrument_recovered';
        riserTag = '.*SPKIR.*\.nc$'; 
        mooringMethod = 'recovered_host';
        riser = load_gc_thredds(mooringSite, riserNode, riserSensor, mooringMethod, riserStream, riserTag);
        riserannotations = get_annotations(mooringSite, riserNode, riserSensor);
        CSPPannotations =  get_annotations(site, node, sensor);
        disp('completed riser load');
        ce = load_gc_thredds(site, node, sensor, method, stream, tag);
        disp('completed profiler load.');
        if nsite==1
            save('CE01ISSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        else
            save('CE07SHSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        end
        disp(int2str(nsite));
    end
end
disp(' Loaded THREDDS');

%% add profile variable to THREDDS
% no CTD data nsite= 2 deployment= 35
% no CTD data nsite= 3 deployment= 6
% no SPKIR data nsite= 3 deployment= 17
% no CTD data nsite= 4 deployment= 15
if 1==addProfile
    cd(folder);
    for nsite = sites
        if nsite==1
            load CE01ISSPspkir.mat;
            ctd=load('CE01ISSP');
        elseif nsite ==2
            load CE02SHSPspkir.mat;
            ctd=load('CE02SHSP');
        elseif nsite==3
            load CE06ISSPspkir.mat;
            ctd=load('CE06ISSP');
        else
            load CE07SHSPspkir.mat;
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
                    legend({'SPKIR','CTD'},'Location','east');
                else
                    disp(['no SPKIR data nsite= ',int2str(nsite),' deployment= ',int2str(i)]);
                end
            else
                disp(['no CTD data nsite= ',int2str(nsite),' deployment= ',int2str(i)]);
            end
        end
        if nsite==1
            save('CE01ISSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        else
            save('CE07SHSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        end
        close all;
    end
end
disp('  Added profile variable');

%ce.Properties.VariableNames' % view variable names

%% look at qc flags

if inspectQC
    cd(folder);
    for nsite = sites
        if nsite==1
            load CE01ISSPspkir.mat;
            ctd=load('CE01ISSP');
        elseif nsite ==2
            load CE02SHSPspkir.mat;
            ctd=load('CE02SHSP');
        elseif nsite==3
            load CE06ISSPspkir.mat;
            ctd=load('CE06ISSP');
        else
            load CE07SHSPspkir.mat;
            ctd=load('CE07SHSP');
        end
        tmp = hour(ce.profiler_datetime);
        ce.day = tmp*0;
        ce.day(tmp <2 | 13 < tmp) = 1; clear tmp;
        % dex=find(ce.spkir_abj_cspp_downwelling_vector<0);
        % ce.spkir_abj_cspp_downwelling_vector(dex)=-ce.spkir_abj_cspp_downwelling_vector(dex);
        for i=1:1:max(ce.deployment)
            figure(double(i));
            dex = find(ce.deployment==i & ~ce.day);
            night_depth = ce.depth(dex);
            night_profiler_datetime = ce.profiler_datetime(dex);
            night_spkir_abj_cspp_downwelling_vector = ce.spkir_abj_cspp_downwelling_vector(dex,:);

            dex = find(ce.deployment==i & ce.day);
            depth = ce.depth(dex);
            profiler_datetime = ce.profiler_datetime(dex);
            spkir_abj_cspp_downwelling_vector = ce.spkir_abj_cspp_downwelling_vector(dex,:);

            if ~isempty(dex)
                subplot(221);
                scatter(profiler_datetime,depth,12,(spkir_abj_cspp_downwelling_vector(:,4)),'filled')
                ylabel('day SPKIR'); hold on; set(gca,'ydir','rev'); 
                a=colorbar; set(a,'limits',[0 1000]);
                % dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                % if ~isempty(dex)
                %     scatter(riser.Time(dex),dex*0+7,10,(riser.spkir_abj_cspp_downwelling_vector(dex,4)));
                % end

                subplot(222);
                dex=find(spkir_abj_cspp_downwelling_vector(:,4)>0);
                scatter(profiler_datetime(dex),depth(dex),12,log10(spkir_abj_cspp_downwelling_vector(dex,4)),'filled')
                ylabel('log 10 day SPKIR no neg'); hold on; set(gca,'ydir','rev'); 
                a=colorbar; set(a,'limits',[0 3]);
                % dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                % if ~isempty(dex)
                %     scatter(riser.Time(dex),dex*0+7,10,(riser.spkir_abj_cspp_downwelling_vector(dex,4)));
                % end

                subplot(223);
                scatter(night_profiler_datetime,night_depth,12,(night_spkir_abj_cspp_downwelling_vector(:,4)),'filled')
                ylabel('night SPKIR'); hold on; set(gca,'ydir','rev'); 
                a=colorbar; set(a,'limits',[0 100]);
                % dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                % if ~isempty(dex)
                %     scatter(riser.Time(dex),dex*0+7,10,(riser.spkir_abj_cspp_downwelling_vector(dex,4)));
                % end
                subplot(224);
                dex=find(night_spkir_abj_cspp_downwelling_vector(:,4)>0);
                scatter(night_profiler_datetime(dex),night_depth(dex),12,log10(night_spkir_abj_cspp_downwelling_vector(dex,4)),'filled')
                ylabel('log10 night SPKIR no neg'); hold on; set(gca,'ydir','rev'); 
                a=colorbar; set(a,'limits',[0 2]);
                % dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                % if ~isempty(dex)
                %     scatter(riser.Time(dex),dex*0+7,10,(riser.spkir_abj_cspp_downwelling_vector(dex,4)));
                % end
            else
                disp(['no SPKIR data in site: ',int2str(nsite),' deployment: ',int2str(i)]);
            end
            set(gcf,'units','inches','position',[1 1 16 10]);
        end

        % if nsite==1
        %     save('CE01ISSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % elseif nsite ==2
        %     save('CE02SHSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % elseif nsite==3
        %     save('CE06ISSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % else
        %     save('CE07SHSPspkir.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % end
    end
end

toc