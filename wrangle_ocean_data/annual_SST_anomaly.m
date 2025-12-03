% get matrix of monthly SST anomaly and annual mean SST by AMM7 gridcell and year

% EJT Dec 2025 tidied up for manuscript

outfpath = 'data';
saveData = 0; 

% read AMM7 mesh
% note these data are available from Copernicus -  see README.txt
filepath = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\AMM7 masks\';
filename_bathy = [filepath 'NWS-MFC_004_001_mask_bathy.nc'];
AMM7.lon = ncread(filename_bathy,'longitude');
AMM7.lat = ncread(filename_bathy,'latitude');
AMM7.deptho = ncread(filename_bathy,'deptho');
[AMM7.LON,AMM7.LAT] = meshgrid(AMM7.lon,AMM7.lat);

n.lon = length(AMM7.lon);
n.lat = length(AMM7.lat);

% monthly temperature output from Copernicus -  see README.txt
datapath = 'D:\NWS_multiyear\phy\thetao'
clear fname

yrs=1993:2022;
n.yrs = length(yrs);
n.months = 12;

% initialise output
annualSST = nan(n.lon,n.lat,n.yrs);
monthlySST = nan(n.lon,n.lat,n.months,n.yrs);

vec_idx = 0;

% loop over years
for i=1:n.yrs
    
    yr=yrs(i)

    YYYY=num2str(yr)

    for mo=1:12
        MM=sprintf('%02d',mo)
        % get monthly AMM7 temperature file
        fname = fullfile(datapath,'monthly',YYYY,['metoffice_foam1_amm7_NWS_TEM_mm' YYYY MM '.nc']);
        temperature = ncread(fname,'thetao');
        monthlySST(:,:,mo,i) = temperature(:,:,1); % surface temperature
    end

    annualSST(:,:,i) = mean( squeeze(monthlySST(:,:,:,i)) ,3);
end



% get monthly anomaly
% gives array of size lon * lat * 12
mo_mean = mean( monthlySST, 4  ); % mean over years for each month; dim 4 is year
mo_anom_SST = monthlySST - mo_mean;

% save to file - need to save mean too so can calculate anomaly in the
% projection data

if saveData
    save(fullfile(outfpath,"annual_AMM7_SST_1993_2022.mat"),"annualSST")
    save(fullfile(outfpath,"monthly_AMM7_SST_anomaly_1993_2022.mat"),"monthlySST","mo_anom_SST","mo_mean");
end