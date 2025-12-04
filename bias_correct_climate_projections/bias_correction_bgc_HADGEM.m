%% get relative change fields
% Calculate monthly mean relative differences (deltas) between bgc
% variables for now and future in the climate projections. deltas will then
% be applied to "now" in the reanalysis to generate a bias-corrected
% (relative change) future time series. Need to calculate: - monthly - for
% near future (2040s) and far future (2090s) - by grid-cell - for all three
% projection datasets - for all predictors

% EJT Dec 2024 updated for nitrate
% EJT Dec 2025 tidied up for manuscript

clear, close all

% bgc projections for HADGEM projection - these files can be downloaded from https://zenodo.org/records/3953801
infpath = 'D:\ECOWINGS\ROAM_HADGEM2\monthly\ersem\';

% where to save deltas
outfpath = 'bias_correct_climate_projections\';

% get mask info - note vertical grids differ by projection
clear M
fname_mesh = 'D:\ECOWINGS\ROAM_HADGEM2\mesh_AMM7.nc';
M.landMask3D = ncread(fname_mesh,'tmask'); % land mask of 0,1
M.vertGridIntervalT = ncread(fname_mesh,'e3t'); % depths of cells

% read lateral T grid from mask file
M.lon_T = ncread(fname_mesh,'glamt');  % lon at T points
M.lat_T = ncread(fname_mesh,'gphit');  % lat at T point
M.lon_T = wrapTo180(M.lon_T); % longitudes need wrapping
n.lons = size(M.lon_T,1);
n.lats = size(M.lat_T,2);

% get water depth
z = M.landMask3D.*M.vertGridIntervalT;
M.depth = cumsum( cat(3,zeros(n.lons,n.lats),z(:,:,1:(end-1))) ,3); % depth at each grid-cell
M.bathymetry=sum(z,3); % sum of grid thicknesses
M.bathymetry(M.bathymetry==0)=NaN;

% make mask of 150 m for NO3 calculations
M.layer_150m=nan(n.lons,n.lats);
for ln=1:n.lons
    for lt=1:n.lats
        if (M.landMask3D(ln,lt,1)~=0)
        idx=find(M.depth(ln,lt,:)<150); % deepest cell within 150 m
        M.layer_150m(ln,lt) = max(idx);
        end
    end
end

% decadal slices
time_period = {'now','NF','FF'};
time_period_yrs = {2010:2023,2040:2049,2090:2099};
n.tp = length(time_period);

% get list of monthly projection files - careful, there are months missing
flist = dir([infpath 'amm7_1m_*ptrc*']); % chl for phytoplankton size fractions are all in here
n.files = length(flist);

% get year of each file
file_year = nan(n.files,1);
for i=1:n.files % better without a loop but works
    file_year(i) = str2double(flist(i).name(9:12));
end

Pfrac = {'chl','diato','dino','pico','nano'};
clear monthly_mean winter_no3

%% get monthly climatology for each time period
% loop over time periods
for k=1:n.tp

    yrs = time_period_yrs{k};
    n.yrs=length(yrs);
    disp(['getting data for ' num2str(yrs(1)) 's'])

    % grab files in date range
    idx=find(file_year>=min(yrs) & file_year<=max(yrs));
    flist_tp=flist(idx); % subset to required years
    n.files_tp = length(flist_tp);

    % get chl  by month for all years (one month per file)
    clear chl chl_all no3 no3_mo
    blank = nan(n.lons,n.lats,n.files_tp);
    for p=1:length(Pfrac)
        chl_all.(Pfrac{p}) = blank;
    end
 
    no3_mo = blank; % monthly 0-150 m no3

    for i=1:n.files_tp
        fname = fullfile(infpath,flist_tp(i).name);
        disp(flist_tp(i).name)

        % --- chl
        chl.diato = ncread(fname,'Chl1');
        chl.dino = ncread(fname,'Chl2');
        chl.pico = ncread(fname,'Chl3');
        chl.nano = ncread(fname,'Chl4');

        % change fill value to NaN
        fillValue = max(chl.diato,[],"all");
        for P = ["diato","dino","nano","pico"]
            chl.(P)(chl.(P)==fillValue)=NaN;
        end

        chl.total=chl.diato+chl.dino+chl.nano+chl.pico;

        % integrate over depth
        chl.total_depthInt = sum( (chl.total .* M.vertGridIntervalT), 3,'omitmissing')./M.bathymetry;
        chl.diato_depthInt = sum( (chl.diato .* M.vertGridIntervalT), 3,'omitmissing')./M.bathymetry;
        chl.dino_depthInt = sum( (chl.dino .* M.vertGridIntervalT), 3,'omitmissing')./M.bathymetry;
        chl.pico_depthInt = sum( (chl.pico .* M.vertGridIntervalT), 3,'omitmissing')./M.bathymetry;
        chl.nano_depthInt = sum( (chl.nano .* M.vertGridIntervalT), 3,'omitmissing')./M.bathymetry;
    
        chl_all.chl(:,:,i) = chl.total_depthInt;

        % calculate deltas using absolute values, not fractions
        chl_all.diato(:,:,i) = chl.diato_depthInt;%./chl.total_depthInt;
        chl_all.dino(:,:,i) = chl.dino_depthInt;%./chl.total_depthInt;
        chl_all.pico(:,:,i) = chl.pico_depthInt;%./chl.total_depthInt;
        chl_all.nano(:,:,i) = chl.nano_depthInt;%./chl.total_depthInt;

        % --- no3

        no3 = ncread(fname,'N3n'); % 297x375x51
        fillValue = ncreadatt(fname,"N3n","_FillValue");
        idx=find(no3==fillValue);
        if ~isempty(idx)
            disp('FOUND FILLERS')
        end

        % integrate over 0-150m
        no3_upper = nan(n.lons,n.lats);
        for ln=1:n.lons
            for lt=1:n.lats
                z150=M.layer_150m(ln,lt);
                if ~isnan(z150)
                    no3_upper(ln,lt) = sum( (no3(ln,lt,1:z150) .* z(ln,lt,1:z150)), 3,'omitmissing')./sum( z(ln,lt,1:z150) );
                end
            end
        end
        no3_mo(:,:,i) = no3_upper;
    end

    % get no3_winter, mean over Jan-Mar
    no3_Jan =  no3_mo(:,:,1:12:n.files_tp);
    no3_Feb =  no3_mo(:,:,2:12:n.files_tp);
    no3_Mar =  no3_mo(:,:,3:12:n.files_tp);
    no3_winter.(time_period{k}) = (no3_Jan+no3_Feb+no3_Mar)/3;

    % get mean chl and no3 for each month
    for p=1:length(Pfrac)
        monthly_mean.(Pfrac{p}).(time_period{k}) = nan(n.lon,n.lat,12);

        for j=1:12
            monthly_mean.(Pfrac{p}).(time_period{k})(:,:,j) = mean(chl_all.(Pfrac{p})(:,:,j:12:n.files_tp),3);
        end
    end

    monthly_mean.no3_mo.(time_period{k}) = nan(n.lon,n.lat,12);
    for j=1:12
        monthly_mean.no3_mo.(time_period{k})(:,:,j) = mean(no3_mo(:,:,j:12:n.files_tp),3);
    end

end

%% (2) get NF/now and FF/now
% (deltas calculated using relative change for bgc variables which must be non-negative)
clear Delta
for k=2:n.tp
    for Pfrac = ["chl","diato","dino","pico","nano"]
        Delta.(Pfrac).(time_period{k}) = monthly_mean.(Pfrac).(time_period{k}) ./ monthly_mean.(Pfrac).now;
    end
end

% monthly no3
for k=2:n.tp
    Delta.no3_mo.(time_period{k}) = monthly_mean.no3_mo.(time_period{k}) ./ monthly_mean.no3_mo.now;
end

% winter no3: take mean over years first
for k=2:n.tp
    Delta.no3_winter.(time_period{k}) = mean(no3_winter.(time_period{k}),3,'omitmissing') ./ mean(no3_winter.now,3,'omitmissing');
end

%% save monthly bgc Deltas
save(fullfile(outfpath,'Delta_CHL_PFT_NO3_HADGEM'),"Delta")
