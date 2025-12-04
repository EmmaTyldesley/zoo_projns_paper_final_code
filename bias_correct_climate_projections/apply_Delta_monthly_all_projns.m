%% apply Delta corrections to AMM7 reanalysis monthly output
% save ready for transforming into R random forest df format

% EJT modified Dec 2024 to include nitrate
% EJT modified Mar 2025 for general tidy up
% EJT modified Dec 2025 for manuscipt

clear, close all, clc
outfpath = 'bias_correct_climate_projections';

% which climate projection model?
% note: same for all except need to load bloom phenology for RECICLE (no
% daily bgc output for HADGEM so cannot calculate bloom metrics)
mdlName = 'HADGEM';%'RECICLE_IPSL'; % 'RECICLE_GFDL'

% --- read AMM7 reanalysis mesh and get depth info needed for integration
% note: this data can be obtained through Copernicus Marine Data:
% https://data.marine.copernicus.eu/product/NWSHELF_MULTIYEAR_BGC_004_011/
% and https://data.marine.copernicus.eu/product/NWSHELF_MULTIYEAR_PHY_004_009
filepath = 'D:\NWS_multiyear\';
filename_bathy = [filepath 'masks\NWS-MFC_004_001_mask_bathy.nc'];
AMM7.lon = ncread(filename_bathy,'longitude');
AMM7.lat = ncread(filename_bathy,'latitude');
AMM7.deptho = ncread(filename_bathy,'deptho');
AMM7.mask = ncread(filename_bathy,'mask');
AMM7.depth = ncread(filename_bathy,'depth');

n.lon = length(AMM7.lon);
n.lat = length(AMM7.lat);
n.z = length(AMM7.depth);

zlayers = AMM7.depth;
zlayers(1) = 1; % surface layer is up to 1 m
dz = diff([0 ; zlayers]); % depth of each band
dz_mat = nan(n.lon,n.lat,n.z);

% create grid cell depths 3D array
for i=1:length(AMM7.lon)
    for j=1:length(AMM7.lat)
        dz_mat(i,j,:) = dz.*squeeze(AMM7.mask(i,j,:));
    end
end
dz_mat_sum=sum(dz_mat,3);


%% load Deltas for climate projection model
% calculated by bias_correction_<variables>_<model_name>.m
% amend to read combined file
infpath = 'bias_correct_climate_projections';

load(fullfile(infpath,['Delta_SST_NBT_SSS_MLD_updated_' mdlName '.mat'])); % Delta.<var>.<NF/FF> lon*lat*12
Delta_all = Delta;
load(fullfile(infpath,['upper_ocean_salinity_updated_' mdlName '.mat'] )); % mean for each time period
load(fullfile(infpath,['Delta_CHL_PFT_NO3_' mdlName '.mat'])); % Delta.<var>.<NF/FF> lon*lat*12

% combine bgc and phys Deltas
P=fieldnames(Delta);
for i=1:length(P)
    Delta_all.(P{i}) = Delta.(P{i});
end

time_periods = {'NF','FF'};
n.tp=length(time_periods);

%% load AMM7 "now" (2010-2022)
% these data can be ontained from
% https://data.marine.copernicus.eu/product/NWSHELF_MULTIYEAR_PHY_004_009/
% and https://data.marine.copernicus.eu/product/NWSHELF_MULTIYEAR_BGC_004_011/
amm7path = 'D:\NWS_multiyear\' %bgc/chl/monthly/';% then <phy/bgc>/<var>/<year>/<file>
%clear predictors
yrs = 2010:2023;
n.yrs = length(yrs);
n.mos = 12;

vars = {'sst','sss','mld','no3_mo','chl'}%,'diato','dino','pico','nano'};
n.vars = length(vars);
varsFileName = {'TEM','SAL','MLD','NITR','CPWC'}%,'DIATO_CPWC','DINO_CPWC','PICO_CPWC','NANO_CPWC'};
varsName = {'thetao','so','mlotst','no3','chl'}%,'diato','dino','pico','nano'};

for i=1:n.vars
    disp(['getting ' vars{i}])
    if any(strcmp(vars{i},{'diato','dino','pico','nano'}))
        searchStr = [amm7path '*/*/monthly/*/*/metoffice_foam1_amm7_NWS_' varsFileName{i} '_mm*.nc']; % mean monthly data
    else
        searchStr = [amm7path '*/*/monthly/*/metoffice_foam1_amm7_NWS_' varsFileName{i} '_mm*.nc']; % mean monthly data

    end
    flist = dir(searchStr);
    n.files = length(flist);

    % get date of file. surely don't need a loop?
    yearStamp = nan(n.files,1);
    for j=1:n.files
        yearStamp(j) = str2double(flist(j).name((end-8):(end-5)));
    end
    % get files in year range
    idx=find(yearStamp>=min(yrs) & yearStamp<=max(yrs));
    flist=flist(idx); % subset to required years
    n.files_in = length(flist);

    % initialise collated variable
    var_all = nan(n.lon,n.lat,n.mos,n.yrs);

    for j=1:n.files_in
        % need to collate into lon*lat*mo*year
        fname = fullfile(flist(j).folder,flist(j).name);
        disp(flist(j).name)
        t = ncread(fname,'time');
        t = datetime(t,'ConvertFrom','epochtime','Epoch','1970-01-01 00:00:00');
        % get year number
        ynum = find(yrs==year(t));
        % get month number
        mnum = month(t);

        % for sst and sss, just read surface
        if strcmp(vars{i},'sst') | strcmp(vars{i},'sss')
            var_t = ncread(fname,varsName{i},[1 1 1 1],[Inf Inf 1 Inf]);
        else
            % for chl and mld, read full depth
            var_t = ncread(fname,varsName{i});
        end
        % for chl, depth integrate
        if any(strcmp(vars{i},{'chl','diato','dino','pico','nano'}))
            test = var_t.*dz_mat;
            test_sum = sum( test, 3,'omitnan');
            chl_depth_int = test_sum./dz_mat_sum;
            var_t = chl_depth_int;
        end

        % for no3, integrate over 0-150 m
        if any(strcmp(vars{i},{'no3_mo'}))
            var_t=integrateByDepth(zlayers,1,150,var_t);
        end

        var_all(:,:,mnum,ynum) = var_t;
    end

    % collate with other vars
    predictors.(vars{i}).now = var_all;
end

%% get winter no3
predictors.no3_winter.now = squeeze( mean( predictors.no3_mo.now(:,:,1:3,:),3 , 'omitmissing') );

%% get SST anom
% (read again as need whole period 1993-2023 to calculate climatology)
disp('getting SST anomaly')

searchStr = [amm7path '*/*/monthly/*/metoffice_foam1_amm7_NWS_TEM_mm*.nc']; % mean monthly temperature
flist = dir(searchStr);
n.files = length(flist);
yrs_full = 1993:2023;

sst_all = nan(n.lon,n.lat,n.mos,length(yrs_full));
for i=1:n.files
    % collate into lon*lat*mo*year
    fname = fullfile(flist(i).folder,flist(i).name);
    disp(flist(i).name)
    t = ncread(fname,'time');
    t = datetime(t,'ConvertFrom','epochtime','Epoch','1970-01-01 00:00:00');
    % get year number
    ynum = find(yrs_full==year(t));
    % get month number
    mnum = month(t);

    sst_t = ncread(fname,'thetao',[1 1 1 1],[Inf Inf 1 Inf]);
    sst_all(:,:,mnum,ynum) = sst_t;
end

sst_mo_mean_1993_2023 = mean(sst_all,4);
sst_anom_all = sst_all - sst_mo_mean_1993_2023; % this looks okay

% grab 2010-2023 as "now"
idx_t = find(ismember(yrs_full,yrs));
predictors.mo_sst_anom.now = sst_anom_all(:,:,:,idx_t);

%% get monthly no3 anomaly wrt 1993-2023
% have to read again as need whole period 1993-2023 to calculate
% climatology
disp('getting no3 anomaly')

searchStr = [amm7path '*/*/monthly/*/metoffice_foam1_amm7_NWS_NITR_mm*.nc']; % mean monthly no3
flist = dir(searchStr);
n.files = length(flist);
yrs_full = 1993:2023;
% initialise collated variable
no3_all = nan(n.lon,n.lat,n.mos,length(yrs_full));

for i=1:n.files
    % collate into lon*lat*mo*year
    fname = fullfile(flist(i).folder,flist(i).name);
    disp(flist(i).name)
    t = ncread(fname,'time');
    t = datetime(t,'ConvertFrom','epochtime','Epoch','1970-01-01 00:00:00');

    ynum = find(yrs_full==year(t));% get year & month number
    mnum = month(t);
    
    if ~isempty(ynum)
        no3_t = ncread(fname,'no3');
        no3_t = integrateByDepth(zlayers,1,150,no3_t);
        no3_all(:,:,mnum,ynum) =no3_t;
    end
end

no3_mo_mean_1993_2023 = mean(no3_all,4); % get monthly climatology
no3_anom_all = no3_all - no3_mo_mean_1993_2023; % get anomaly

% grab 2010-2023 as "now"
idx_t = find(ismember(yrs_full,yrs));
predictors.no3_mo_anom.now = no3_anom_all(:,:,:,idx_t);

%% apply monthly Deltas to make adjusted time series
blank = nan(n.lon,n.lat,n.mos,n.yrs);

for i=1:n.vars % loop over variables
    for j=1:n.tp % loop over time periods
        % initialise adjusted variable
        var_future = blank;
        for k=1:n.yrs
            if any(strcmp(vars{i},{'sst','nbt'})) % additivie
                var_future(:,:,:,k) = predictors.(vars{i}).now(:,:,:,k) + Delta_all.(vars{i}).(time_periods{j})(:,:,:);
            else % multiplicative
                var_future(:,:,:,k) = predictors.(vars{i}).now(:,:,:,k) .* Delta_all.(vars{i}).(time_periods{j})(:,:,:);  
            end
        end
        % bung into structure
        predictors.(vars{i}).(time_periods{j}) = var_future;
    end
end

%% apply monthly Deltas to no3_winter
    for j=1:n.tp % loop over time periods

var_future = nan(n.lon,n.lat,n.yrs);
        for k=1:n.yrs
                var_future(:,:,k) = predictors.no3_winter.now(:,:,k) .* Delta_all.no3_winter.(time_periods{j});

        end
        % bung into structure
        predictors.no3_winter.(time_periods{j}) = var_future;
    end

%% convert phytoplankton chl components to fractions
% only if have included them
if any(strcmp('diato',vars))

    Pfrac={'diato','dino','pico','nano'}
    for i=1:n.tp
        for p=1:length(Pfrac)
            predictors.(Pfrac{p}).(time_periods{i}) = predictors.(Pfrac{p}).(time_periods{i}) ./ predictors.chl.(time_periods{i});
        end
    end
end

%% calculate mo_sst_anom wrt 1993-2023
for i=1:n.tp
    predictors.mo_sst_anom.(time_periods{i}) = predictors.sst.(time_periods{i}) - sst_mo_mean_1993_2023;
end

%% calculate no3_mo_anom wrt 1993-2023
for i=1:n.tp
    predictors.no3_mo_anom.(time_periods{i}) = predictors.no3_mo.(time_periods{i}) - no3_mo_mean_1993_2023;
end

%% do the same for the annual SPG proxy
% load historical SPG proxy
SPG_proxy_past = readtable('C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\ECOWINGS\data\SPG_proxy_annual.txt');

% get mean 1993-2023 for calculating anomaly
mean_SPG_proxy_past = mean(SPG_proxy_past.SPG_proxy_annual);

% get 2010-2022 section
SPG_proxy_now = SPG_proxy_past.SPG_proxy_annual(SPG_proxy_past.yr>=min(yrs) & SPG_proxy_past.yr<=max(yrs))
SPG_proxy_anom_now = SPG_proxy_now - mean_SPG_proxy_past;

% apply Delta
for i=1:n.tp
    Delta_all.SPG_proxy_annual.(time_periods{i})= upper_ocean_salinity.(time_periods{i}) - upper_ocean_salinity.now;
    SPG_proxy_future.(time_periods{i}) = SPG_proxy_now + Delta_all.SPG_proxy_annual.(time_periods{i});
    SPG_proxy_anom_future.(time_periods{i}) =  SPG_proxy_future.(time_periods{i}) - mean_SPG_proxy_past;
end

figure
tiledlayout(2,1)
nexttile
plot(yrs,SPG_proxy_now)
hold on
for i=1:n.tp
    plot(yrs,SPG_proxy_future.(time_periods{i}))
end
grid on
title('annual upper ocean salinity')
plot(SPG_proxy_past.yr,SPG_proxy_past.SPG_proxy_annual,'--')
legend('now','NF','FF','full past','location','southwest')

nexttile
plot(yrs,-SPG_proxy_anom_now)
hold on
for i=1:n.tp
    plot(yrs,-SPG_proxy_anom_future.(time_periods{i}))
end
grid on

title('SPG proxy: annual upper ocean salinity anomaly wrt 1993-2018')
print(fullfile(outfpath,'SPG_proxy_projns_HADGEM.pdf'),'-bestfit','-dpdf')

% combine with other 
predictors.SPG_proxy_anom = SPG_proxy_anom_future;
predictors.SPG_proxy_anom.now = SPG_proxy_anom_now;

%% load bias corrected bloom projections
if ~strcmp(mdlName,'HADGEM')
    temp=load(fullfile(infpath,'output',['bias_corrected_' mdlName '_predictors_bloom.mat']));
    predictors.bloom_start = temp.predictors.bloom_start;
    predictors.bloom_duration = temp.predictors.bloom_duration;
end
%% save
fname = fullfile(outfpath,['bias_corrected_' mdlName '_predictors_incl_no3.mat']);
disp(['saving output: ' fname])
save(fname,"predictors",'-v7.3') % large file