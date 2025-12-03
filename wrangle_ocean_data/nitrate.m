%% calculate metrics of nitrate from amm7 reanalysis: 
% annual mean upper ocean pre-bloom (Jan-Mar)
% and monthly upper ocean nitrate anomaly (wrt 1993-2023 as for SST)
% cf. Skagseth 2022, Hatun 2017

% EJT tidied Dec 2025 for manuscript submission

clear, close all

% load monthly amm7 nitrate
infpath = 'D:\NWS_multiyear\bgc\no3\monthly';
% note the above data are available from https://data.marine.copernicus.eu/product/NWSHELF_MULTIYEAR_PHY_004_009/description
% example filename:
fname = fullfile(infpath,'2012','metoffice_foam1_amm7_NWS_NITR_mm201201.nc');

outfpath = 'data';

yrs = 1993:2023;
n.yrs = length(yrs);

zlayers = ncread(fname,'depth');
ztop=0; 
zbot=150; % integrate 0-150 m

%% get upper ocean no3 by year, month and gridcell
no3_upper_monthly=nan(297,375,n.yrs,12);

for i=1:n.yrs
    YYYY = num2str(yrs(i));
    disp(YYYY)
    for mo=1:12 % Jan to Mar

        MM = sprintf('%02d',mo);
        disp(MM)

        fname = fullfile(infpath,YYYY, ...
            ['metoffice_foam1_amm7_NWS_NITR_mm' YYYY MM '.nc']);

        no3=ncread(fname,'no3');
        no3_upper_monthly(:,:,i,mo) = integrateByDepth(zlayers,ztop,zbot,no3);

    end
end

% pre-bloom no3 by year and gridcell
no3_upper_winter = mean(no3_upper_monthly(:,:,:,1:3),4,'omitmissing'); % mean over Jan-Mar each year

%% get monthly anomaly
% gives array of size lon * lat * 12
mo_mean = mean( no3_upper_monthly, 3  ); % mean over years for each month; dim 3 is year
no3_mo_anom = no3_upper_monthly - mo_mean;

%% save
no3_yrs=yrs;
dataname = 'amm7_monthly_nitrate_1993_2023.mat';
save(fullfile(outfpath,dataname),"no3_upper_winter","no3_upper_monthly","no3_mo_anom","no3_yrs")
