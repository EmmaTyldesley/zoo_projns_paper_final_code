%% get relative change fields - RECICLE
% Calculate monthly mean differences (deltas) between now and future in the
% climate projections physical output. Deltas will then be added to "now"
% in the reanalysis to generate a bias-corrected (relative change) future
% time series. Need to calculate: - monthly - for near future (2040s) and
% far future (2090s) - by grid-cell - for all three projection datasets -
% for all predictors

% updated Nov 2024 to include DEB driver near bottom temperature (NBT)
% updated Feb 2025 to make NBT and SST additive, SSS and MLD relative

% EJT Dec 2025 - updated for manuscript

clear, close all
doCalc=1;

% --- get RECICLE mesh
% note this can be obtained from BODC at
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

% --- get SPG proxy location (highest corr with SPG idx)
lon0=-12.1112;
lat0=54.0674;
[I0,J0]=findnearest(lon0,lat0,M.lon_T(:,1),M.lat_T(1,:));
% get depth bands - will need this for integration
z = squeeze(M.depth(I0,J0,:)); % depth layers at yur, yol
idx_upper_ocean = find(z>50 & z<500);   % bands within 50-500 m
dz = diff(z(idx_upper_ocean));

% --- filepaths to monthly data
% note that daily output can be obtained from BODC links given above &
% converted to monthly using NCO's ncra command or similar (see README.txt)
mdlName = 'RECICLE_IPSL';% 'RECICLE_IPSL'; %'RECICLE_GFDL';

if strcmp(mdlName,'RECICLE_IPSL')
    infpath = 'D:\RECICLE_IPSL-CM5A-MR\bodc\PML230161\';
else
    infpath = 'D:\RECICLE_GFDL-ESM2G\bodc\PML230162\';
end
outfpath = 'bias_correct_climate_projections\';

% --- which time periods
time_period = {'now','NF','FF'};
time_period_yrs = {2010:2022,2040:2049,2090:2099};
n.tp = length(time_period);

clear monthly_mean
vars = {'sst','nbt','sss','mld'}; % need to add chl code in here
n.vars=length(vars);

% critical T change for MLD calculation
deltaT_crit = 0.2; % AMM7 hindcast uses 0.2

% there is an error in file AMM7_1m_201105_grid_T_from_daily.nc
% time stamp for 5/5/2011 is 1/1/2050 and all data missing
% load the daily file and recalculate
if strcmp(mdlName,'RECICLE_IPSL')
    fname = 'D:\RECICLE_IPSL-CM5A-MR\bodc\PML230161\1990-2039\2011\AMM7_1d_201105_grid_T.nc';
    salinity_depth = squeeze( ncread(fname,'Sal',[I0 J0 1 1],[1 1 Inf Inf]) );
    salinity_depth = mean(salinity_depth(:,[1:4 6:end]),2); % remove the zero column and take mean over remaining months
    SPG_proxy_201105 = sum( salinity_depth(idx_upper_ocean(1:(end-1))).*dz ) / sum(dz);
    S_201105 = ncread(fname,'Sal');
    S_201105 = mean(S_201105(:,:,:,[1:4 6:end]),4);
    T_201105 = ncread(fname,'Temp');
    T_201105 = mean(T_201105(:,:,:,[1:4 6:end]),4);
end

% get list of monthly files
close all
figure
if doCalc
    flist = dir([infpath '*\*\*_grid_T_from_daily*']); % use the generated daily-monthly files
    n.files = length(flist);
    % get year of each file
    file_year = nan(n.files,1);
    file_mo = nan(n.files,1);
    for i=1:n.files
        file_year(i) = str2double(flist(i).name(9:12));
        file_mo(i) = str2double(flist(i).name(13:14));
    end

    % -- make mask of bottom layer
    % read in example T
    fname =fullfile(flist(1).folder,flist(1).name);
    T = ncread(fname,'Temp');
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
        file_year_tp=file_year(idx);
        file_mo_tp=file_mo(idx);
        n.files_tp = length(flist_tp);

        % -- get variables by month for all years (one month per file)

        % initialise
        clear sst sss nbt mld upper_ocean_sal no3_winter no3_mo
        blank = nan(297,375,n.files_tp);
        sst = blank; % sea surface temperature
        nbt = blank; % near bottom temperature
        sss = blank; % sea surface salinity
        mld = blank; % mixed layer depth

        upper_ocean_sal=nan(n.files_tp,1);

        for i=1:n.files_tp
            fname = fullfile(flist_tp(i).folder,flist_tp(i).name);
            disp(flist_tp(i).name)

            % need to read full depth T and S to calc MLD and NBT
            T = ncread(fname,'Temp');
            S = ncread(fname,'Sal');
            if strcmp(mdlName,'RECICLE_IPSL') & file_year_tp(i)==2011 & file_mo_tp(i)==5
                disp(['correcting T and S for 2011 05'])
                T = T_201105;
                S = S_201105;
            end

            sst(:,:,i) = T(:,:,1);
            sss(:,:,i) = S(:,:,1);

            for ln=1:n.lons
                for lt=1:n.lats
                    if ~isnan(bottom_layer(ln,lt))
                        nbt(ln,lt,i) = T(ln,lt,bottom_layer(ln,lt));
                    end
                end
            end

            % read upper ocean salinity
            upper_ocean_sal_depth = squeeze( ncread(fname,'Sal',[I0 J0 1 1],[1 1 Inf Inf],[1 1 1 1]) ); % read surface layer only
            upper_ocean_sal_depth=upper_ocean_sal_depth(idx_upper_ocean(1:(end-1)));
            upper_ocean_sal(i) = sum( upper_ocean_sal_depth.*dz ) / sum(dz);
            if strcmp(mdlName,'RECICLE_IPSL') & file_year_tp(i)==2011 & file_mo_tp(i)==5
                disp(['correcting SPG proxy for 2011 05'])
                upper_ocean_sal(i)=SPG_proxy_201105;
            end


            % now calculate density
            % note using this to calculate mixed layer depth using
            % definition consistent with the reanalysis model
            [rho,rho_diff]=seawater_density(S,T,0);
            rho = rho_diff; % just call it rho for clarity (i.e. density - 1000)
            clear rho_diff

            T_0 = sst(:,:,i); % reference T
            S_0 = sss(:,:,i); % reference S and rho at same depth
            rho_0 = squeeze(rho(:,:,1));

            T_z = T_0-deltaT_crit;

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
                            %I=I(I>1); % don't allow it to pick the surface layer
                            mld(ln,lt,i) = M.depth(ln,lt,I(1));
                        else
                            mld(ln,lt,i) = M.depth(ln,lt,end); % mixed water column
                        end
                    end
                end
            end

        end

        % take mean for each month
        % don't do this for SPG proxy - using annual
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

    % (2) get NF-now and FF-now : sst and nbt additive; mld and sss relative
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
    save(fullfile(outfpath,['Delta_SST_NBT_SSS_MLD_updated_' mdlName]),"Delta")
    save(fullfile(outfpath,['upper_ocean_salinity_updated_' mdlName]),"upper_ocean_salinity")
    save(fullfile(outfpath,['Delta_additive_SSS_' mdlName]),"Delta_add_SSS") % save this too - useful information

end