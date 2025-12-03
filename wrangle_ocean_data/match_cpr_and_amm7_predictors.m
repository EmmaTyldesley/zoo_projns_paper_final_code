%% build table of cpr samples and predictors for use in randomForest
% loops over cpr samples and matches with nearest AMM7 grid-cell in space
% and time predictors: x, y, SST, annual SST, monthly SST anomaly, SSS, mld, yearday, chl,
% diato_frac, dino_frac, day/night, water depth, annual winter NO3, monthly NO3.
% Other predictors added
% later by R code: SPG index, Norw. Sea RFC

% EJT updated Jan 2024 to include mld & rerun on updated cpr dataset
% EJT updated Dec 2024 to include nitrate
% EJT tidied Dec 2025 for manuscript submission

clear, clc, close all
cd("C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\ECOWINGS\zoo_projns_paper_final_code\");
outfpath = 'output';

saveData=0; % not saving, currently

%% load cpr zooplankton energy (kJ/m3)
load('data\cpr_add.mat') % this includes additional species sent by D Johns
cpr=cpr_add;
cpr.yday = day(cpr.t,'dayofyear');
cpr.yr = year(cpr.t);
n.samples = length(cpr.x);

%% read AMM7 reanalysis mesh
filepath = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\AMM7 masks\';
filename_bathy = [filepath 'NWS-MFC_004_001_mask_bathy.nc'];
% note, this mask file is available from
% https://data.marine.copernicus.eu/product/NWSHELF_MULTIYEAR_PHY_004_009/description
AMM7.lon = ncread(filename_bathy,'longitude');
AMM7.lat = ncread(filename_bathy,'latitude');
AMM7.deptho = ncread(filename_bathy,'deptho');

%% load previously collated AMM7 reanalysis data
load('data\annual_AMM7_SST_1993_2022.mat'); % 
annualSST_yrs = 1993:2022;
load('data\monthly_AMM7_SST_anomaly_1993_2022.mat'); % monthly SST anomaly 1993-2022
load('data\amm7_monthly_nitrate_1993_2023.mat'); % upper (0-150m) ocean nitrate - monthly, monthly anomaly, and annual mean pre-winter (Jan-Mar)
load('data\bloom_1998_2022.mat') % phytoplankton bloom dynamics

%% loop over cpr samples, finding nearest AMM7 grid-cell
% slow - only run if cpr2amm7_idx.mat missing from data folder
blank = nan(n.samples,1);
if 0
    cpr2amm7.lon=blank;
    cpr2amm7.lat=blank;
    %samples_in_AMM7 = find( cpr.x>=-20.0 & cpr.x<=13.0 & cpr.y>=40.0 & cpr.y<=65.0);
    for i=1:n.samples
        i
        % only for samples within amm7 domain
        if cpr.x(i)>=-20.0 & cpr.x(i)<=13.0 & cpr.y(i)>=40.0 & cpr.y(i)<=65
            [cpr2amm7.lon(i),cpr2amm7.lat(i)]=findnearest(cpr.x(i),cpr.y(i),AMM7.lon,AMM7.lat);
        end
    end
    % save
    save(fullfile("data","cpr2amm7_idx.mat"),'cpr2amm7')
else
    load(fullfile("data","cpr2amm7_idx.mat"))
end

%% get depth
amm7_predictors.depth = blank;
for i=1:n.samples
    % only for samples within amm7 domain
    if ~isnan(cpr2amm7.lon(i)) & ~isnan(cpr2amm7.lat(i)) 
        amm7_predictors.depth(i) = AMM7.deptho(cpr2amm7.lon(i),cpr2amm7.lat(i)); 
    end
end

%% get AMM7 predictors

% initialise predictors
predictors = {'sst','annual_sst','mo_sst_anom','sss','chl_amm7','mld','diato','dino','nano','pico','bloom_start','bloom_duration','no3_mo_anom','no3_winter'};
for i=1:length(predictors)
    amm7_predictors.(predictors{i})=blank;
end

% set depth bands to nan
zlayers = nan;

% loop over samples and grab AMM7 values
fdir = 'D:\NWS_multiyear\'; % note the data are available from https://data.marine.copernicus.eu/product/NWSHELF_MULTIYEAR_PHY_004_009/description
for i=1:n.samples

    yr = year(cpr.t(i));
    mo = month(cpr.t(i));

    if yr>=1993 & ~isnan(cpr2amm7.lon(i)) & ~isnan(cpr2amm7.lat(i)) % check if in AMM7 reanalysis domain

        YYYY = num2str(yr);
        disp(YYYY)
        MM = sprintf('%02d',month(cpr.t(i)));
        DD = sprintf('%02d',day(cpr.t(i)));

        start = [cpr2amm7.lon(i) cpr2amm7.lat(i) 1 1];   % location in nc variable
        count = [1 1 1 1]; % just want that one value
        count_depth = [1 1 Inf 1]; % for getting full depth (chl)

        % sst
        fname.sst = fullfile(fdir,'phy','thetao/daily',YYYY,MM,['metoffice_foam1_amm7_NWS_TEM_dm' YYYY MM DD '.nc']);
        amm7_predictors.sst(i) = ncread(fname.sst,'thetao',start,count); % just read the value we want, not whole variable

        % annual SST
        idx_yr = find(yr==annualSST_yrs);
        amm7_predictors.annual_sst(i) = annualSST(cpr2amm7.lon(i),cpr2amm7.lat(i),idx_yr);

        % monthly SST anomaly
        amm7_predictors.mo_sst_anom(i) = mo_anom_SST(cpr2amm7.lon(i),cpr2amm7.lat(i),mo,idx_yr);

        % salinity
        fname.sss = fullfile(fdir,'phy','so/daily',YYYY,MM,['metoffice_foam1_amm7_NWS_SAL_dm' YYYY MM DD '.nc']);
        amm7_predictors.sss(i) = ncread(fname.sss,'so',start,count); 

        % chl - depth integrated
        fname.chl = fullfile(fdir,'bgc','chl/daily',YYYY,MM,['metoffice_foam1_amm7_NWS_CPWC_dm' YYYY MM DD '.nc']);

        % read depth bands if not got already
        if isnan(zlayers)
            zlayers = ncread(fname.chl,'depth');
            zlayers(1) = 1; % surface layer is up to 1 m
            dz = diff([0 ; zlayers]); % height of each depth band
        end

        chl_with_depth = squeeze( ncread(fname.chl,'chl',start,count_depth) ); % read full depth & do weighted sum

        % integrate over depth
        chl_depth_int = sum( chl_with_depth .* dz , 'omitnan') / sum(dz(~isnan(chl_with_depth)));
        amm7_predictors.chl_amm7(i) = chl_depth_int;
         
        % mixed layer depth
        fname.mld = fullfile(fdir,'phy','mld/daily',YYYY,MM,['metoffice_foam1_amm7_NWS_MLD_dm' YYYY MM DD '.nc']);
        amm7_predictors.mld(i)  = ncread(fname.mld,'mlotst',[cpr2amm7.lon(i) cpr2amm7.lat(i) 1],[1 1 1]); 

        % phytoplankton functional types
        % this is just surface
        for pft = ["DIATO","DINO","NANO","PICO"]
            fname.pft = fullfile(fdir,'bgc/pft/daily',pft,YYYY,MM,'metoffice_foam1_amm7_NWS_'+pft+'_CPWC_dm'+YYYY+MM+DD+'.nc');
            amm7_predictors.(lower(pft))(i) = ncread(fname.pft,lower(pft),start,count); 
            %amm7_predictors.(lower(pft))(i) = P(cpr2amm7.lon(i),cpr2amm7.lat(i),1,1); % 'mg chl m^{-3}'
        end

        % bloom phenology
        % only for 1998 on
        if yr>=1998
            idx = find(bloom.years==yr);
            amm7_predictors.bloom_start(i) = bloom.start(cpr2amm7.lon(i),cpr2amm7.lat(i),idx);
            amm7_predictors.bloom_duration(i) = bloom.duration(cpr2amm7.lon(i),cpr2amm7.lat(i),idx);
        end

        % nitrate
        idx=find(no3_yrs==yr);
        amm7_predictors.no3_winter(i) = no3_upper_winter(cpr2amm7.lon(i),cpr2amm7.lat(i),idx); % annual pre-bloom
        amm7_predictors.no3_mo_anom(i) = no3_mo_anom(cpr2amm7.lon(i),cpr2amm7.lat(i),idx,mo); % monthly

    end
end

% save struct
% if saveData
%     save(fullfile(outfpath,'amm7_predictors.mat'),'amm7_predictors');
% end

% turn into table and save for reading into R
zed_and_amm7 = [struct2table(cpr) struct2table(amm7_predictors)];

% convert diato, dino, nano and pico conc into frac of total
P_total = zed_and_amm7.diato+zed_and_amm7.dino+zed_and_amm7.nano+zed_and_amm7.pico;
zed_and_amm7.diato_frac = zed_and_amm7.diato ./ P_total;
zed_and_amm7.dino_frac = zed_and_amm7.dino ./ P_total;
zed_and_amm7.nano_frac = zed_and_amm7.nano ./ P_total;
zed_and_amm7.pico_frac = zed_and_amm7.pico ./ P_total;

if saveData
    writetable(zed_and_amm7, fullfile('outputs','zed_and_amm7.txt'));
end

%% check plots - this all looks good
close all
% check geographical variation
figure(1), clf
scatter(cpr.x,cpr.y,6,log10(amm7_predictors.no3_winter))
% check seasonality
figure(2), clf
plot(cpr.t,amm7_predictors.no3_mo_anom,'o')



