%% get relative change fields
% Calculate monthly mean differences (deltas) between now and future in the
% climate projections physical output. Deltas will then be added to "now"
% in the reanalysis to generate a bias-corrected (relative change) future
% time series. Need to calculate: - monthly - for near future (2040s) and
% far future (2090s) - by grid-cell - for all three projection datasets -
% for all predictors

% EJT Dec 2025 - updated for manuscript

clear, close all
doCalc=1;

% get HADGEM mesh note this can be obtained from
% https://zenodo.org/records/3953801
clear M
fname_mesh = 'D:\ECOWINGS\ROAM_HADGEM2\mesh_AMM7.nc';
M.lon_T = ncread(fname_mesh,'glamt');  % lon at T points
M.lat_T = ncread(fname_mesh,'gphit');  % lat at T point
M.lon_T = wrapTo180(M.lon_T); % longitudes need wrapping
n.lons = size(M.lon_T,1);
n.lats = size(M.lat_T,2);

M.landMask3D = ncread(fname_mesh,'tmask'); % land mask of 0,1
M.vertGridIntervalT = ncread(fname_mesh,'e3t'); % depths of cells
z = M.landMask3D.*M.vertGridIntervalT;
M.depth = cumsum( cat(3,zeros(n.lons,n.lats),z(:,:,1:(end-1))) ,3); % depth at each grid-cell
M.bathymetry=sum(z,3); % sum of grid thicknesses
M.bathymetry(M.bathymetry==0)=NaN;

% get salinity-based sub-polar gyre (SPG) proxy location (grid-cell with
% highest correlation with SPG idx) Note: this was established by
% calculating correlation coefficients between the SPG index of
% Hatun-Chafik 2018 (https://doi.org/10.1029/2018JC014101) and upper ocean
% salinity in the Copernicus AMM7 model reanalysis
% (https://data.marine.copernicus.eu/product/NWSHELF_MULTIYEAR_PHY_004_009/description)
lon0=-12.1112;
lat0=54.0674;
[I0,J0]=findnearest(lon0,lat0,M.lon_T(:,1),M.lat_T(1,:));

% get depth bands - need this for integration
z = squeeze(M.depth(I0,J0,:)); % depth layers at location
idx_upper_ocean = find(z>50 & z<500);   % bands within 50-500 m
dz = diff(z(idx_upper_ocean));

% monthly projection data note: can obtain this from Zenodo link given
% above
infpath = 'D:\ECOWINGS\ROAM_HADGEM2\monthly\nemo\';
outfpath = 'bias_correct_climate_projections\'

% (1) read  projection "now" (2010-2022) and get monthly means by gridcell
% ==> array n.lat*n.lon*12 months

time_period = {'now','NF','FF'};
time_period_yrs = {2010:2022,2040:2049,2090:2099};
n.tp = length(time_period);

clear monthly_mean
vars = {'sst','nbt','sss','mld'};
n.vars=length(vars);

% critical T change for mixed layer depth (MLD) calculation
deltaT_crit = 0.2; % AMM7 hindcast uses 0.2. see user manual: https://documentation.marine.copernicus.eu/PUM/CMEMS-NWS-PUM-004-009-011.pdf


%% get list of monthly files
if doCalc
    flist = dir([infpath 'amm7_1m_*_grid_T*']);
    n.files = length(flist);
    % get year of each file
    file_year = nan(n.files,1);
    for i=1:n.files
        file_year(i) = str2double(flist(i).name(9:12));
    end

    % -- make mask of bottom layer read in example T
    fname =fullfile(flist(1).folder,flist(1).name);
    T = ncread(fname,'votemper');
    bottom_layer=nan(n.lons,n.lats);
    for ln=1:n.lons
        for lt=1:n.lats
            above_seabed = find(~isnan(T(ln,lt,:)));
            if ~isempty(above_seabed)
                bottom_layer(ln,lt) = max(above_seabed);
            end
        end
    end

    % loop over time periods
    for k=1:n.tp

        yrs = time_period_yrs{k};
        n.yrs=length(yrs);
        disp(['getting data for ' num2str(yrs(1)) 's'])

        % grab files in date range
        idx=find(file_year>=min(yrs) & file_year<=max(yrs));
        flist_tp=flist(idx); % subset to required years
        n.files_tp = length(flist_tp);

        % get variables by month for all years (one month per file)
        clear sst sss
        blank = nan(297,375,n.files_tp);
        sst = blank; % sea surface temperature
        nbt = blank;
        sss = blank; % sea surface salinity
        mld = blank; % mixed layer depth
        upper_ocean_sal=nan(n.files_tp,1);

        for i=1:n.files_tp
            fname = fullfile(infpath,flist_tp(i).name);
            disp(flist_tp(i).name)
            sst(:,:,i) = ncread(fname,'votemper',[1 1 1 1],[Inf Inf 1 Inf],[1 1 1 1]); % read surface layer only
            sss(:,:,i) = ncread(fname,'vosaline',[1 1 1 1],[Inf Inf 1 Inf],[1 1 1 1]); % read surface layer only

            for ln=1:n.lons
                for lt=1:n.lats
                    if ~isnan(bottom_layer(ln,lt))
                        nbt(ln,lt,i) = T(ln,lt,bottom_layer(ln,lt));
                    end
                end
            end

            % read upper ocean salinity
            upper_ocean_sal_depth = squeeze( ncread(fname,'vosaline',[I0 J0 1 1],[1 1 Inf Inf],[1 1 1 1]) ); % read surface layer only
            upper_ocean_sal_depth=upper_ocean_sal_depth(idx_upper_ocean(1:(end-1)));
            upper_ocean_sal(i) = sum( upper_ocean_sal_depth.*dz ) / sum(dz);

            % also need to read full depth T and S to calc MLD
            T = ncread(fname,'votemper');
            S = ncread(fname,'vosaline');
            % now calculate density
            [rho,rho_diff]=seawater_density(S,T,0);
            rho = rho_diff; % just call it rho for clarity (i.e. density - 1000)
            clear rho_diff

            T_0 = sst(:,:,i); % reference T
            S_0 = sss(:,:,i); % reference S and rho at same depth
            rho_0 = squeeze(rho(:,:,1));

            T_z = T_0-deltaT_crit; % temperature 0.8 deg C lower than ref

            [rho_Tz,rho_diff_Tz] = seawater_density(S_0,T_z,0); % density if T changed by 0.2C but S and p fixed
            rho_Tz = rho_diff_Tz;
            clear rho_diff_Tz

            delta_rho_crit_T = abs(rho_0 - rho_Tz); % critical density; lon*lat*t

            % now loop over lon, lat to get MLD at each grid cell
            for ln=1:n.lons
                for lt = 1:n.lats
                    if ~isnan(M.bathymetry(ln,lt))
                        delta_rho = (  squeeze(rho(ln,lt,:)) - rho_0(ln,lt));
                        I = find ( delta_rho  > delta_rho_crit_T(ln,lt) );
                        if ~isempty(I)
                            mld(ln,lt,i) = M.depth(ln,lt,I(1));
                        else
                            mld(ln,lt,i) = M.depth(ln,lt,end); % mixed water column
                        end
                        if mld(ln,lt,i)==0
                            disp('found one!:')
                        end
                    end
                end
            end

        end

        % now take mean for each month don't do this for SPG proxy - using
        % annual
        for v=1:n.vars
            monthly_mean.(vars{v}).(time_period{k}) = nan(297,375,12);

        end

        for j=1:12
            monthly_mean.sst.(time_period{k})(:,:,j) = mean(sst(:,:,j:12:n.files_tp),3);
            monthly_mean.sss.(time_period{k})(:,:,j) = mean(sss(:,:,j:12:n.files_tp),3);
            monthly_mean.mld.(time_period{k})(:,:,j) = mean(mld(:,:,j:12:n.files_tp),3);
            monthly_mean.nbt.(time_period{k})(:,:,j) = mean(nbt(:,:,j:12:n.files_tp),3);
        end

        % get mean SPG_proxy for time period

        upper_ocean_salinity.(time_period{k}) = mean(upper_ocean_sal);

    end

    % (2) get NF-now and FF-now : sst and nbt additive, sss both, and mld
    % relative
    clear Delta Delta_add_SSS
    for k=2:n.tp % only for NF and FF
        for v=1:n.vars
            if strcmp(vars{v},'sst') | strcmp(vars{v},'nbt') % additive
                Delta.(vars{v}).(time_period{k}) = monthly_mean.(vars{v}).(time_period{k}) - monthly_mean.(vars{v}).now;
            else % relative
                Delta.(vars{v}).(time_period{k}) = monthly_mean.(vars{v}).(time_period{k}) ./ monthly_mean.(vars{v}).now;
            end
        end
        Delta_add_SSS.(time_period{k}) = monthly_mean.sss.(time_period{k}) - monthly_mean.sss.now;

    end

    % save monthly Deltas

    save(fullfile(outfpath,'Delta_SST_NBT_SSS_MLD_updated_HADGEM'),"Delta")
    save(fullfile(outfpath,'upper_ocean_salinity_updated_HADGEM'),"upper_ocean_salinity")
    save(fullfile(outfpath,'Delta_additive_SSS_HADGEM'),"Delta_add_SSS")

end



