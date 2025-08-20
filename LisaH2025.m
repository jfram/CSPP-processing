% 2025-08-20
% OOI Data for Lisa Hildebrand MMI student
clear; close all; clc;
folder = 'C:\Users\jfram\OneDrive - Oregon State University\Documents\MATLAB\CSPPproc';
cd(folder);
load('CE01ISSPdostaT.mat');

% for units, see https://github.com/oceanobservatories/ooi-data-explorations/blob/master/python/examples/notebooks/dosta/ooiea_dosta_quality_assessments.ipynb

dCast=readtable("Endurance_Small_Boat_Discrete_Summary.xlsx","Sheet","Endurance_Small_Boat_Discrete_S");
dCast.depth=dCast.CTDDepth_m_;
dCast.temperature=dCast.CTDTemperature1_degC_;
dCast.conductivity=dCast.CTDConductivity1_S_m_;
dCast.salinity=dCast.CTDSalinity1_psu_;
dCast.dissolved_oxygen=dCast.CTDOxygen_mL_L_; % 32 g/mol 
for i=1:length(dCast.depth)
    tmp=char(dCast.CTDBottleClosureTime_UTC_(i));
    tmp=tmp(1:(end-1));
    dCast.dtime(i)=datetime(tmp,'timezone','utc');
end
dex=find(dCast.Cruise=="EK20250804");
dCast.dissolved_oxygen(dex)=[3.533,3.063,2.714,6.115,5.789,5.657,7.483,7.476,7.490]*44.6596; % from ml/L to umol/L
% readtable wasn't perfect, so I copied the numbers directly out of the spreadsheet


%% all together now
figure
subplot(211);
plot(riser.Time,riser.dosta_ln_optode_oxygen,'.',bottom.Time,bottom.dosta_ln_optode_oxygen,'.',dCast.dtime(dex),dCast.dissolved_oxygen(dex),'gx','MarkerSize',16)
ct=datetime(2025,8,[4 6],'timezone','utc'); ylabel('Dissolved Oxygen (umol/L)'); legend({'NSIF','MFN','discrete'});
xlim(ct); title('NSIF daytime peaks may indicate biofouling: mid-summer. All UTC');
subplot(212);
plot(riser.Time,riser.dissolved_oxygen,'.',bottom.Time,bottom.dissolved_oxygen,'.',ce.Time,ce.dissolved_oxygen,'.','MarkerSize',16)
pt=datetime(2025,6,[27 29],'timezone','utc'); legend({'NSIF','MFN','profiler'});
xlim(pt); title('Profiler data line up with moored data (mostly): early summer');
ylabel('Dissolved Oxygen (umol/L)')

% discrete cast times
dr=find(ct(1) <= riser.Time & riser.Time <= ct(2));
db=find(ct(1) <= bottom.Time & bottom.Time <= ct(2));
dc=find(ct(1) <= dCast.dtime(dex) & dCast.dtime(dex) <= ct(2));
figure
% plot(riser.Time(dr),riser.dissolved_oxygen(dr),'o',bottom.Time(db),bottom.dissolved_oxygen(db),'d',dCast.dtime(dex(dc)),dCast.dissolved_oxygen(dex(dc)),'x','MarkerSize',12,'linewidth',2);
scatter(riser.Time(dr),7+dr*0,30,riser.dosta_ln_optode_oxygen(dr),'filled'); hold on;
scatter(bottom.Time(db),23+db*0,30,bottom.dosta_ln_optode_oxygen(db),'filled');
scatter(dCast.dtime(dex(dc)),dCast.depth(dex(dc)),70,dCast.dissolved_oxygen(dex(dc)),'filled');
set(gca,'ydir','reverse'); colorbar; ylabel('depth (m)'); caxis([0 500]); 
set(gca,'fontsize',14); title('moored and discrete samples: mid-summer'); 

% profile times
dr=find(pt(1) <= riser.Time & riser.Time <= pt(2));
db=find(pt(1) <= bottom.Time & bottom.Time <= pt(2));
dc=find(pt(1) <= ce.Time & ce.Time <= pt(2));
figure
% plot(riser.Time(dr),riser.dissolved_oxygen(dr),'o',bottom.Time(db),bottom.dissolved_oxygen(db),'d',ce.Time(dc),ce.dissolved_oxygen(dc),'.','MarkerSize',12,'linewidth',2);
scatter(riser.Time(dr),7+dr*0,30,riser.dosta_ln_optode_oxygen(dr),'filled'); hold on;
scatter(bottom.Time(db),23+db*0,30,bottom.dosta_ln_optode_oxygen(db),'filled');
scatter(ce.Time(dc),ce.depth(dc),20,ce.estimated_oxygen_concentration(dc),'filled');
set(gca,'ydir','reverse'); colorbar; caxis([0 500]); 
set(gca,'fontsize',14); ylabel('depth (m)');  title('moored and profiler data: early summer'); 