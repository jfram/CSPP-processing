% setup the defaults used to request the data
site = 'CE02SHBP';
node = 'LJ01D';
sensor = '06-CTDBPN106';
method = 'streamed';
stream = 'ctdbp_no_sample';
parallel = true;  % invokes use of parallel processing to handle downloads

% load the common file processing utilities
ph = process_files;

% determine the number of deployments of this particular sensor
deployments = list_deployments(site, node, sensor);

% resampling settings
dt = minutes(15);
vars = {'deployment', 'depth', 'sea_water_pressure', 'sea_water_pressure_qartod_results', ...
    'sea_water_electrical_conductivity', 'sea_water_electrical_conductivity_qartod_results', ...
    'sea_water_temperature', 'sea_water_temperature_qartod_results', ...
    'sea_water_practical_salinity', 'sea_water_practical_salinity_qartod_results', 'sea_water_density', ...
    'corrected_dissolved_oxygen', 'corrected_dissolved_oxygen_qartod_results'};

% for each deployment, download the data and resample to a 15-minute, median
% averaged data record (clean-up some variables, etc. as we go)
data = cell(numel(deployments), 1);  % pre-allocate array to save results in
for i = 1:numel(deployments)
    tag = sprintf('deployment%04d.*CTDBP.*\\.nc$', i);
    ctd = load_gc_thredds(site, node, sensor, method, stream, tag, parallel);
    if isempty(ctd) == false
        ctd = ctd(:, vars);  % sub-select only the variables we want to keep
        m = ctd.sea_water_electrical_conductivity_qartod_results == 4;
        ctd(m==1, :) = []; clear m  % remove out-of-water data
        ctd.Time = ctd.Time + (dt / 2);  % adjust the time to center the retime window
        data{i} = retime(ctd, 'regular', 'median', 'TimeStep', dt);  % append the data
    end %if
    %%% could put the steps below here to save the data per deployment if
    %%% you wanted, or just run everything and have one final dataset
end %for
clear deployments dt vars i tag ctd

% finally, concatenate the timetables together returning the results
if length(data) > 1
    data = ph.merge_frames(data);
else
    data = data{1};
end %if
clear site node sensor method stream parallel ph
