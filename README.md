# CSPP-processing
Process data from the NSF OOI Endurance Array Coastal Surface Piercing Profilers
%% Profile analysis plan
% download CTD data from OOI's THREDDS catalog -- doner
% organize by deployment -- done. in the data stream
% organize by profile -- done (one deployment is missing: CE06ISSP 00006. It is in CI's queue to be fixed)
% plot by profile -- done 
%   note automated QC issues. Not available as of 2024-12-23. They are in CI's queue. Check back early 2025.
%   note annotations -- use get_annotations.m from Wingard
%   make mask for bad data (go through each deployment)
% download NSIF, MFN, and surface CT data -- done
%   decimated CE02 BEP CTD & vel data from Brandy -- done
% check data against moored CTD data -- 
%   plot 221: 1: scatter T at the depth for each deployment. -- done 
%             2: 1:1 T profiler vs. moored (3 dots per profile, so this is per deployment). 
%             3 & 4: are S
%   plot 131: 1: T vs. D for each profile in each deployment. done
%             2: S vs. D. -- done
%             3: T-S profile and moored CTD with depth as color for each profile 
% CTD lags
%   ProcessRealTime?
%   C a few cm down. 12 cm from the top to the center of the C cell.
%   D instant but below. 56 cm from top to the pressure sensor. 
%   T instant and exposed. 1 cm from the top to the T sensor.
%   make variable named depth_sensor for each sensor (depth_temperature, depth_conductivity, depth_pressure). We want value as a function of depth. 
% generate additional annotations (make notes per deployment, insert into CI)
% concatenate deployments into grids: 5cm bin (CTD only), 25cm bin (all sensors), and 2m bin (NH-10 for CTD and ADCP).
% get data from other CSPP instruments.
% depth_sensor 
%   measure average for all profilfers.  
%   compare to QCTs and deployment sheets
%   for DOSTA, see from deployment sheets and QCTs when depth changes from center of the profilers to the top. (Iridium ones by Mari)
% QC other data streams
% check other CSPP data against NSIF and surface data
% correct instruments
%    night time PAR and SPKIR should be zero
%    recalculate nitrate with shifted CT data.
% grid into 25 cm bins (i.e. the requirements)
