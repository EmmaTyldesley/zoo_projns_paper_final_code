%% get relative change fields - bgc - RECICLE
% (phytoplankton bloom phenology calculated separately because daily rather than monthly output required)

% EJT Dec 2024 - modified to include nitrate
% EJT Dec 2025 - updated for manuscript

clear, close all

% ---- which projection
mdlName = 'RECICLE_GFDL'% 'RECICLE_IPSL';

% --- get RECICLE mesh
clear M
% mesh and mask - same for both RECICLE runs
% Note these can be obtained from BODC at
% https://doi.org/10.5285/07877700-0e22-5fb5-e063-6c86abc03058 and
% https://doi.org/10.5285/0786d770-ba54-0e57-e063-6c86abc09fdf
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
M.depth = cumsum(z,3); % depth at each grid-cell
M.bathymetry = squeeze(M.depth(:,:,end));

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

% --- filepath to monthly projection data
% note that daily output can be obtained from BODC links given above &
% converted to monthly using NCO's ncra command or similar (see README.txt)
if strcmp(mdlName,'RECICLE_IPSL')
    infpath = 'D:\RECICLE_IPSL-CM5A-MR\bodc\PML230161/';
else
    infpath = 'D:\RECICLE_GFDL-ESM2G\bodc\PML230162\';
end
outfpath = 'bias_correct_climate_projections';

% --- which time periods
time_period = {'now','NF','FF'};
time_period_yrs = {2010:2022,2040:2049,2090:2099};
n.tp = length(time_period);

% --- get list of monthly files
flist = dir([infpath '*\*\*_bgc_T_from_daily*']); % use the generated monthly files
n.files = length(flist);
% get year of each file
file_year = nan(n.files,1);
for i=1:n.files 
    file_year(i) = str2double(flist(i).name(9:12));
end

% phytoplankton size fractions
Pfrac = {'chl','diato','dino','pico','nano'};

clear monthly_mean no3_winter

%% loop over time periods
for k=1:n.tp

    yrs = time_period_yrs{k};
    n.yrs=length(yrs);
    disp(['getting data for ' num2str(yrs(1)) 's'])

    % grab files in date range
    idx=find(file_year>=min(yrs) & file_year<=max(yrs));
    flist_tp=flist(idx); % subset to required years
    n.files_tp = length(flist_tp);

    % get bgc  by month for all years (one month per file)
    clear chl chl_all no3 no3_mo
    blank = nan(297,375,n.files_tp);

    for p=1:length(Pfrac)
        chl_all.(Pfrac{p}) = blank;
    end

    no3_mo = blank; % monthly 0-150m no3

    for i=1:n.files_tp

        fname = fullfile(flist_tp(i).folder,flist_tp(i).name);
        disp(flist_tp(i).name)

        % --- chl
        % read each chl fraction, sum and depth integrate
        chl.diato = ncread(fname,'P1_Chl'); % diatom chl, mg C.m^{-3}
        chl.dino = ncread(fname,'P4_Chl'); % dinoflagellates
        chl.pico = ncread(fname,'P3_Chl'); % picoP
        chl.nano = ncread(fname,'P2_Chl'); % nanoP

        %fillValue - no fill value for RECICLE
        chl.total=chl.diato+chl.dino+chl.nano+chl.pico;

        % integrate over depth
        chl.total_depthInt = sum( (chl.total .* z), 3,'omitmissing')./M.bathymetry;
        chl.diato_depthInt = sum( (chl.diato .* z), 3,'omitmissing')./M.bathymetry;
        chl.dino_depthInt = sum( (chl.dino .* z), 3,'omitmissing')./M.bathymetry;
        chl.pico_depthInt = sum( (chl.pico .* z), 3,'omitmissing')./M.bathymetry;
        chl.nano_depthInt = sum( (chl.nano .* z), 3,'omitmissing')./M.bathymetry;
    
        chl_all.chl(:,:,i) = chl.total_depthInt;
        % calculate deltas using absolute values, not fractions
        chl_all.diato(:,:,i) = chl.diato_depthInt;%./chl.total_depthInt;
        chl_all.dino(:,:,i) = chl.dino_depthInt;%./chl.total_depthInt;
        chl_all.pico(:,:,i) = chl.pico_depthInt;%./chl.total_depthInt;
        chl_all.nano(:,:,i) = chl.nano_depthInt;%./chl.total_depthInt;

        % --- no3
        no3 = ncread(fname,'N3_n'); % 297x375x51
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

    % now take mean for each month - careful, might need to match months
    for p=1:length(Pfrac)
        monthly_mean.(Pfrac{p}).(time_period{k}) = nan(297,375,12);

        for j=1:12
            monthly_mean.(Pfrac{p}).(time_period{k})(:,:,j) = mean(chl_all.(Pfrac{p})(:,:,j:12:n.files_tp),3);
        end
    end

    monthly_mean.no3_mo.(time_period{k}) = nan(297,375,12);
    for j=1:12
        monthly_mean.no3_mo.(time_period{k})(:,:,j) = mean(no3_mo(:,:,j:12:n.files_tp),3);
    end
end

% (2) get NF/now and FF/now
% use relative change 
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

%% save monthly Deltas
save(fullfile(outfpath,['Delta_CHL_PFT_NO3_' mdlName]),"Delta")
