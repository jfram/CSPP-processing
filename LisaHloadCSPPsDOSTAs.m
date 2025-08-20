clear; close all; clc;
tic
% One last time, load telemetered data

% https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/catalog.html
downloadTHREDDS=1;
addProfile=1;
inspectQC=1;
folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder);
nsitedepths=25;

%% gather THREDDS CSPP data
if downloadTHREDDS
    % set up required inputs -- CE01ISSP
    site = 'CE01ISSP';
    sensor = '02-DOSTAJ000';
    mooringSite = 'CE01ISSM';
    riserNode = 'RID16';
    bottomNode = 'MFD37';
    riserSensor =  '03-DOSTAD000';
    bottomSensor = '03-DOSTAD000';
    riserStream =  'dosta_abcdjm_ctdbp_dcl_instrument';
    bottomStream = riserStream;
    node = 'SP001';
    method = 'recovered_cspp';
    stream = 'dosta_abcdjm_cspp_instrument_recovered';
    tag = 'deployment0025.*DOSTA.*\.nc$';
    riserTag = 'deployment0022.*DOSTA.*\.nc$'; bottomTag = riserTag;
    mooringMethod = 'telemetered';
    riser = load_gc_thredds(mooringSite, riserNode, riserSensor, mooringMethod, riserStream, riserTag);
    riserannotations = get_annotations(mooringSite, riserNode, riserSensor);
    CSPPannotations =  get_annotations(site, node, sensor);

    bottom = load_gc_thredds(mooringSite, bottomNode, bottomSensor, mooringMethod, bottomStream, bottomTag);
    bottomannotations = get_annotations(mooringSite, bottomNode,bottomSensor);
    ce = load_gc_thredds(site, node, sensor, method, stream, tag);

    save('CE01ISSPdostaT.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
end
disp(' Loaded THREDDS');

% ce.Properties.VariableNames' % view variable names

%% add profile variable to THREDDS
if 1==addProfile
    cd(folder);
    load CE01ISSPdostaT.mat;
    ctd=load('CE01ISSP');
    ce.profiler_datetime = datetime(ce.profiler_timestamp,'ConvertFrom','posixtime');
    ce.profiler_datetime.TimeZone = 'UTC';
    ce.profile = ce.deployment*NaN;
    % find profile number
    for i=25
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
                disp(['no DOSTA data nsite= ',int2str(1),' deployment= ',int2str(i)]);
            end
        else
            disp(['no CTD data nsite= ',int2str(1),' deployment= ',int2str(i)]);
        end
    end
    save('CE01ISSPdostaT.mat','-v7.3','ce','riser','bottom','riserannotations','bottomannotations','CSPPannotations');
end
disp('  Added profile variable');

toc