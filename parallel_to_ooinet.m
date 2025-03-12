% setup the defaults used to request the data
site = 'CE02SHBP';
node = 'LJ01D';
sensor = '06-CTDBPN106';
method = 'streamed';
stream = 'ctdbp_no_sample';
dataset = strjoin({site, node, sensor, method, stream}, "/");
ooinet_defaults
dap_url = "http://opendap.oceanobservatories.org/thredds/fileServer";
ph = process_files;

% determine the number of deployments of this particular sensor
deployments = list_deployments(site, node, sensor);

% resampling settings
dt = minutes(15);
vars = {'corrected_dissolved_oxygen', 'corrected_dissolved_oxygen_qartod_results', 'deployment', ...
    'depth', 'sea_water_electrical_conductivity', 'sea_water_electrical_conductivity_qartod_results', ...
    'sea_water_temperature', 'sea_water_temperature_qartod_results', 'sea_water_pressure', ...
    'sea_water_pressure_qartod_results', 'sea_water_practical_salinity', ...
    'sea_water_practical_salinity_qartod_results', 'sea_water_density'};

% create a parallel processor pool using the defaults based on the machine
pobj = parpool('local', 5);
data = cell(max(deployments), 1);  % pre-allocate an array to hold the final data set

% request the data
for d = 1:numel(deployments)
    % get the start and stop dates for each deployment
    [strt, stop] = get_deployment_dates(site, node, sensor, deployments(d));
    
    % create a date range, with a 30 day step
    date_start = strt:30:stop-30;
    date_stop = strt+30:30:stop;
    if strt == 0 && stop == 0
        continue
    end %if
    
    % use a parallel for loop to index through the date range, requesting and
    % downloading 30 day chunks of the data
    deploy = cell(size(date_start, 2));
    parfor i = 1:size(date_start, 2)
        dstrt = datetime(date_start(i), 'ConvertFrom', 'datenum', 'Format',"yyyy-MM-dd'T'HH:mm:ss.000'Z'");
        dstop = datetime(date_stop(i), 'ConvertFrom', 'datenum', 'Format',"yyyy-MM-dd'T'HH:mm:ss.000'Z'");
        nc_files = M2M_Call(string(dataset), string(dstrt), string(dstop), options);
        if numel(nc_files) >= 1
            temp = cell(numel(nc_files), 1)
            for j = 1:numel(nc_files)
                nc_url = join([dap_url, '/', nc_files(j), '#fillmismatch'], '')
                source = websave([tempname, '.nc'], nc_url, weboptions("Timeout", 300));
                ctd = ph.process_file(source); %#ok<PFBNS>
                if isempty(ctd) == false
                    ctd = ctd(:, vars);  % sub-select only the variables we want to keep
                    m = ctd.sea_water_electrical_conductivity_qartod_results == 4;
                    ctd(m==1, :) = [];  % remove out-of-water data
                    ctd.Time = ctd.Time + (dt / 2);  % adjust the time to center the retime window
                    temp{j} = retime(ctd, 'regular', 'median', 'TimeStep', dt);  % append the data
                end %if
                delete(source);
            end %for
            if length(temp) > 1
                temp = ph.merge_frames(temp);
            else
                temp = temp{1};
            end %if
        end %if
        deploy{i} = temp;
    end %parfor
    data{d} = ph.merge_frames(deploy);
end %for
clear d strt stop date_start date_stop i dstrt dstop nc_files temp j nc_url ctd m

% shutdown the parallel processor pool
delete(pobj)
clear pobj

% wrap it all up 
data = ph.merge_frames(data);
