%% --- get bias corrected projection bloom phenology ---
% calculate annual bloom phenology metrics from climate projections
% calculate Delta change values
% apply to amm7 reanalysis
% note: can only do for RECICLE projections - no daily resolution HADGEM BGC

% updated EJT March 2025 for:
% - bug (calculating relative Delta then wrongly
% applying additive)
% - general tidy up
% - consistent output file path

% updated EJT Dec 2025 for manuscript

clear, clc, close all

doCalc=0; % takes a while so save and reload once calculated

% have to run this one model at a time because need harddrive connected
mdlName = 'RECICLE_IPSL'; %'RECICLE_GFDL';

% --- get RECICLE mesh
% note these can be obtained from BODC at
% https://doi.org/10.5285/07877700-0e22-5fb5-e063-6c86abc03058 and
% https://doi.org/10.5285/0786d770-ba54-0e57-e063-6c86abc09fdf
clear M
% mesh and mask info - same for both RECICLE runs
fname_mesh_h = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\ECOWINGS\RECICLE mesh files\mesh_hgr.nc';
fname_mask = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\ECOWINGS\RECICLE mesh files\mask.nc';
fname_mesh_z = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\ECOWINGS\RECICLE mesh files\mesh_zgr.nc';
% read T grid. note lon doesn't need wrapping.
M.lon_T = ncread(fname_mesh_h,'glamt');  % lon at T points
M.lat_T = ncread(fname_mesh_h,'gphit');  % lat at T point
n.lons = size(M.lon_T,1);
n.lats = size(M.lat_T,2);
% read bathymetry and depth layers
M.landMask3D = ncread(fname_mask,'tmask'); % land mask of 0,1
M.vertGridIntervalT = ncread(fname_mesh_z,'e3t_0'); % has time varying surface, hence _0
z = double(M.landMask3D).*M.vertGridIntervalT;
M.depth = cumsum( cat(3,zeros(n.lons,n.lats),z(:,:,1:(end-1))) ,3);
fname_bathy = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\ECOWINGS\RECICLE mesh files\bathy_meter.nc';
M.bathymetry = ncread(fname_bathy,'Bathymetry');

% --- filepaths to daily bgc
% obtained from BODC. see links above.
if strcmp(mdlName,'RECICLE_IPSL')
    infpath = 'D:\RECICLE_IPSL-CM5A-MR\bodc\PML230161\';
else
    infpath = 'D:\RECICLE_GFDL-ESM2G\bodc\PML230162\';
end
outfpath = 'bias_correct_climate_projections\';
plotpath='bias_correct_climate_projections\'

% get list of **daily** files
flist = dir([infpath '/*/*/*_bgc_T.nc']); % e.g. AMM7_1d_201504_bgc_T.nc
n.files = length(flist);

time_periods = {'now','NF','FF'};
n.tp=3;
time_periods_years = {2010:2023,2040:2049,2090:2099};

%% --- Step 1: calculate bloom metrics
% bloom start = day of 15% cumulative chl
% bloom duration = days between 15% and 85% cumulative chl

if doCalc % only run if not already calculated
    for p=1:n.tp % loop over time periods
        disp(['getting bloom dynamics for ' time_periods{p}])

        yrs = time_periods_years{p};
        n.yrs = length(yrs);

        %initialise bloom metrics
        blank = nan(297,375,length(yrs)); % lon, lat, year
        bloom_start.(time_periods{p}) = blank;
        bloom_duration.(time_periods{p}) = blank;
        bloom_stop.(time_periods{p}) = blank;

        % loop over years and grab chl
        for i=1:n.yrs
            % initalise chl array
            chl_annual = nan(n.lons,n.lats,366);
            rowcount=1;

            % data are stored in two chunks
            if yrs(i)<2040
                insubfolder = '1990-2039';
            else
                insubfolder = '2040-2099'
            end

            % loop over months
            for mo=1:12

                % load file, integrate chl with depth, concatenate
                fname = fullfile(infpath,insubfolder,num2str(yrs(i)),['AMM7_1d_' num2str(yrs(i)) sprintf('%02d',mo) '_bgc_T.nc']);
                disp(fname)

                % read chl fractions and sum
                % this is slow
                chl.diato = ncread(fname,'P1_Chl'); % diatom chl, mg C.m^{-3}
                chl.dino = ncread(fname,'P4_Chl'); % dinoflagellates
                chl.pico = ncread(fname,'P3_Chl'); % picoP
                chl.nano = ncread(fname,'P2_Chl'); % nanoP
                % fillValue - no fill value for RECICLE
                chl.total=chl.diato+chl.dino+chl.nano+chl.pico;

                % integrate over depth
                chl.total_depthInt = sum( (chl.total .* M.vertGridIntervalT), 3,'omitmissing')./M.bathymetry;

                % concatenate
                n.days = eomday(yrs(i),mo);
                chl_annual(:,:,rowcount:rowcount+(n.days-1)) = squeeze(chl.total_depthInt);
                rowcount = rowcount+n.days;
            end

            % calculate bloom metrics for this year
            disp('got all the data. now calculating bloom dynamics...')
            daysinyear = sum( eomday( yrs(i), 1:12 ) );
            t_annual = 1:daysinyear;
            chl_cumsum = cumsum(chl_annual(:,:,1:daysinyear),3);      % get cumulative sum along t dimension
            start_conc = 0.15*chl_cumsum(:,:,end);    % P conc for start of bloom
            end_conc = 0.85*chl_cumsum(:,:,end);      % P conc for end of bloom

            for j=1:n.lons
                for k=1:n.lats
                    idx = find( chl_cumsum(j,k,:)>=start_conc(j,k) & chl_cumsum(j,k,:)<=end_conc(j,k) );
                    if ~isempty(idx)
                        bloom_start.(time_periods{p})(j,k,i) = t_annual(idx(1));
                        bloom_stop.(time_periods{p})(j,k,i) = t_annual(idx(end));
                        bloom_duration.(time_periods{p})(j,k,i) = length(idx); % days
                    end
                end
            end
        end
    end
    % save to file
    save(fullfile(outfpath,['bloom_phenology_' mdlName '.mat']),"bloom_start","bloom_stop","bloom_duration");

else
    % or load if already calculated
    load(fullfile(outfpath,['bloom_phenology_' mdlName '.mat']));
end

%% --- Step 2: get multiplicative Delta change values
clear Delta Delta_add

% get mean phenology metrics by time period
for i=1:n.tp
    mean_bloom_start.(time_periods{i}) = mean( bloom_start.(time_periods{i}) ,3,'omitmissing');
    mean_bloom_duration.(time_periods{i}) = mean( bloom_duration.(time_periods{i}) ,3,'omitmissing');
end

% get relative Delta, to avoid negative values
for i=2:n.tp 
    Delta.bloom_start.(time_periods{i}) = mean_bloom_start.(time_periods{i}) ./ mean_bloom_start.now;
    Delta.bloom_duration.(time_periods{i}) = mean_bloom_duration.(time_periods{i}) ./ mean_bloom_duration.now;
end

% get additive Delta too as might need
for i=2:n.tp 
    Delta_add.bloom_start.(time_periods{i}) = mean_bloom_start.(time_periods{i}) - mean_bloom_start.now;
    Delta_add.bloom_duration.(time_periods{i}) = mean_bloom_duration.(time_periods{i}) - mean_bloom_duration.now;
end

% save output
save(fullfile(outfpath,['Delta_bloom_updated_' mdlName '.mat']),"Delta","Delta_add");

%% Step 3: apply Deltas to AMM7 reanalysis 2010-2019
bloom_reanalysis = load(['data\bloom_1998_2022.mat'],'bloom');
idt=find(bloom_reanalysis.bloom.years>=2010 & bloom_reanalysis.bloom.years<=2019);

% get 2010-2019 time series by grid-cell
predictors.bloom_start.now = bloom_reanalysis.bloom.start(:,:,idt);
predictors.bloom_duration.now = bloom_reanalysis.bloom.duration(:,:,idt);

% multiply by Delta values to get future time series
for i=2:n.tp

    predictors.bloom_start.(time_periods{i}) = predictors.bloom_start.now .* Delta.bloom_start.(time_periods{i});
    predictors.bloom_duration.(time_periods{i}) = predictors.bloom_duration.now .* Delta.bloom_duration.(time_periods{i});

    % set any negative bloom dates to 1 (very few, is an edge effect)
    idx=find(predictors.bloom_start.(time_periods{i})<0);
    if(~isempty(idx))
        disp('warning: found negative bloom dynamics')
        predictors.bloom_start.(time_periods{i})(idx)=1;
    end
end

% save bias corrected bloom metrics
fname = fullfile(outfpath,['bias_corrected_' mdlName '_predictors_bloom.mat']);
disp(['saving output: ' fname])
save(fname,"predictors")

%% make test plots
close all
figure
tiledlayout(3,1)
for i=1:n.tp
    nexttile
    histogram(predictors.bloom_start.(time_periods{i}))
    xlim([1 200])
    title(time_periods{i})
    if i==n.tp
        xlabel('yearday')
    end
    grid on
end

sgtitle('Bloom start, whole NE Atl - bias corrected RECICLE IPSL')
%print(fullfile(outfpath,'bloom_start_hist_IPSL.pdf'),'-bestfit','-dpdf')


%% make example time series plot at Dogger Bank & Wee Bankie & somewhere offshore
siteNames = {'Stonehaven','Dogger Bank','Rockall Trough'}
lon0=[-2.113 2.2189  -15]
lat0=[56.9635 54.8575  54.5]

figure
tiledlayout(3,2)
for j=1:length(siteNames)
    [I0,J0]=findnearest(lon0(j),lat0(j),M.lon_T(:,1),M.lat_T(1,:));
    nexttile

    for i=1:n.tp
        plot(squeeze(predictors.bloom_start.(time_periods{i})(I0,J0,:)))
        hold on
    end
    title(['bloom start ' siteNames{j}])

    xlim([1 10])
    xlabel('year of decade')
    ylabel('yearday')
    nexttile
    for i=1:n.tp
        plot(squeeze(predictors.bloom_duration.(time_periods{i})(I0,J0,:)))
        hold on
    end
    if j==3
        legend(time_periods,'location','northeast')
    end
    title(['bloom duration ' siteNames{j}])
    xlim([1 10])
    xlabel('year of decade')

    ylabel('days')
end
sgtitle('Example projected RECICLE GFDL bloom phenology')
%print(fullfile(outfpath,'bloom_examples_GFDL.pdf'),'-bestfit','-dpdf')
