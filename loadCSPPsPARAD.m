clear; close all; clc;
tic
% https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/catalog.html
downloadTHREDDS=0;
addProfile=1;
inspectQC=1;
folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder);
nsitedepths=[25,80,29,87];
sites=1; %1:4

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
        riser = load_gc_thredds(mooringSite, riserNode, riserSensor, mooringMethod, riserStream, riserTag);
        riserannotations = get_annotations(mooringSite, riserNode, riserSensor);
        CSPPannotations =  get_annotations(site, node, sensor);
        disp('completed riser load');
        ce = load_gc_thredds(site, node, sensor, method, stream, tag);
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

%% add profile variable to THREDDS
% no PARAD data nsite= 2 deployment= 27
% no CTD data nsite= 2 deployment= 35
% no CTD data nsite= 3 deployment= 6
% no CTD data nsite= 4 deployment= 15
if 1==addProfile
    cd(folder);
    for nsite = sites
        if nsite==1
            load CE01ISSPparad.mat;
            ctd=load('CE01ISSP');
        elseif nsite ==2
            load CE02SHSPparad.mat;
            ctd=load('CE02SHSP');
        elseif nsite==3
            load CE06ISSPparad.mat;
            ctd=load('CE06ISSP');
        else
            load CE07SHSPparad.mat;
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
                    legend({'PARAD','CTD'},'Location','east');
                else
                    disp(['no PARAD data nsite= ',int2str(nsite),' deployment= ',int2str(i)]);
                end
            else
                disp(['no CTD data nsite= ',int2str(nsite),' deployment= ',int2str(i)]);
            end
        end
        if nsite==1
            save('CE01ISSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite ==2
            save('CE02SHSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        elseif nsite==3
            save('CE06ISSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        else
            save('CE07SHSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
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
            load CE01ISSPparad.mat;
            ctd=load('CE01ISSP');
        elseif nsite ==2
            load CE02SHSPparad.mat;
            ctd=load('CE02SHSP');
        elseif nsite==3
            load CE06ISSPparad.mat;
            ctd=load('CE06ISSP');
        else
            load CE07SHSPparad.mat;
            ctd=load('CE07SHSP');
        end
        tmp = hour(ce.profiler_datetime);
        ce.day = tmp*0;
        ce.day(tmp <2 | 13 < tmp) = 1; clear tmp;
        for i=1:1:max(ce.deployment)
            figure(double(i));
            dex = find(ce.deployment==i & ~ce.day);
            night_depth = ce.depth(dex);
            night_profiler_datetime = ce.profiler_datetime(dex);
            night_parad_j_par_counts_output = ce.parad_j_par_counts_output(dex);
            night_parad_j_par_counts_output_qartod_results = ce.parad_j_par_counts_output_qartod_results(dex);

            dex = find(ce.deployment==i & ce.day);
            depth = ce.depth(dex);
            profiler_datetime = ce.profiler_datetime(dex);
            parad_j_par_counts_output = ce.parad_j_par_counts_output(dex);
            parad_j_par_counts_output_qartod_results = ce.parad_j_par_counts_output_qartod_results(dex);

            if ~isempty(dex)
                subplot(211);
                scatter(profiler_datetime,depth,12,log10(parad_j_par_counts_output),'filled')
                ylabel('log10(PAR)'); hold on; set(gca,'ydir','rev'); 
                a=colorbar; %set(a,'limits',[-1 3.2]);
                % dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
                % if ~isempty(dex)
                %     scatter(riser.Time(dex),dex*0+7,10,riser.spkir_abj_cspp_downwelling_vector(dex,4));
                % end
                dex=find(parad_j_par_counts_output_qartod_results>1);
                if ~isempty(dex)
                    plot(profiler_datetime(dex),parad_j_par_counts_output_qartod_results(dex),'r.','MarkerSize',12);
                end
                title('CSPP log10(PAR) with red QARTOD');

                subplot(212);
                scatter(night_profiler_datetime,night_depth,12,log10(night_parad_j_par_counts_output),'filled')
                ylabel('log10(PAR)'); hold on; set(gca,'ydir','rev'); 
                a=colorbar; %set(a,'limits',[-1 3.2]);
                dex=find(night_parad_j_par_counts_output_qartod_results>1);
                if ~isempty(dex)
                    plot(night_profiler_datetime(dex),night_parad_j_par_counts_output_qartod_results(dex),'r.','MarkerSize',12);
                end
                title('Night CSPP PAR with red QARTOD');            
            else
                disp(['no data in site: ',int2str(nsite),' deployment: ',int2str(i)]);
            end
            % 
            % subplot(212);
            % dex = find(ce.deployment==i);
            % if ~isempty(dex)
            %     dex=find(min(profiler_datetime)<riser.Time & riser.Time<max(profiler_datetime));
            %     plot(riser.Time(dex),riser.dissolved_oxygen(dex),'b-'); hold on;
            %     mintime=datetime(year(min(profiler_datetime)),month(min(profiler_datetime)),day(min(profiler_datetime)));
            %     maxtime=datetime(year(max(profiler_datetime)),month(max(profiler_datetime)),day(max(profiler_datetime)));
            %     dt=.25; % days
            %     tmptime=mintime:dt:(dt+maxtime);  tmptime.TimeZone = 'UTC';
            %     tmpoxygen=(1:length(tmptime))*NaN;
            %     for j=1:length(tmptime)
            %         ind=find(tmptime(j)<=profiler_datetime & profiler_datetime<(tmptime(j)+dt)); % & depth < 5 & depth < 10);
            %         if ~isempty(ind)
            %             tmpoxygen(j)=mean(dissolved_oxygen(ind));
            %         end
            %     end
            %     plot(tmptime,tmpoxygen,'go-'); grid on;
            %     if ~isempty(dex)
            %         legend({'NSIF','CSPP'});
            %     end
            %     title('DO time series NSIF depth');
            %     ylim([0 400]);
            % end

            set(gcf,'units','inches','position',[1 1 16 10]);
        end

        % if nsite==1
        %     save('CE01ISSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % elseif nsite ==2
        %     save('CE02SHSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % elseif nsite==3
        %     save('CE06ISSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % else
        %     save('CE07SHSPparad.mat','-v7.3','ce','riser','riserannotations','CSPPannotations');
        % end
    end
end

toc