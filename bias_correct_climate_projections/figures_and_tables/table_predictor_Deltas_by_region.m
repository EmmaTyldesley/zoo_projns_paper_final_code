%% make table of projection monthly Deltas for each predictor by region

clear, close all, clc
mdlName = {'HADGEM','RECICLE_IPSL','RECICLE_GFDL'};
mdlName_long = {'HADGEM','IPSL','GFDL'};
n.mdls=length(mdlName);

plotpath = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/RF zooplankton climate change paper/plots/';
infpath = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/MATLAB/ECOWINGS/climate projections bias correction/output';

% load Delta values - everything except bloom phenology
clear Delta_all Delta_add Delta_bloom_add % the latter is for additive bloom phenology which is more intuitive to report than % change
load(fullfile(infpath,'Delta_all_models_updated.mat'))

% load bloom phenology Deltas
% note these are annual so lon*lat
% other Deltas are by month so lon*lat*mo
for m=1:n.mdls
    if ~(strcmp(mdlName{m},'HADGEM'))
        load(fullfile(infpath,['Delta_bloom_updated_' mdlName{m} '.mat']));
        Delta_bloom.(mdlName{m})= Delta;
        Delta_all.(mdlName{m}).bloom_start = Delta.bloom_start;
        Delta_all.(mdlName{m}).bloom_duration = Delta.bloom_duration;

        Delta_bloom_add.(mdlName{m}).bloom_start = Delta_add.bloom_start;
        Delta_bloom_add.(mdlName{m}).bloom_duration = Delta_add.bloom_duration;
    end
end

clear Delta Delta_bloom

vars = fieldnames(Delta_all.RECICLE_GFDL);
vars_long = {'sea surface temperature','near bottom temperature','sea surface salinity', ...
    'mixed layer depth','chlorophyll','','','','','nitrate anomaly','winter nitrate','bloom start','bloom duration'};
vars_units = {[char(0176) 'C'],[char(0176) 'C'],'','log_{10} m','mg.m^{-3}', ...
    '','','','','mmol.m^{-3}','mmol.m^{-3}','yearday','days'};
vars_delta_units = {[char(0176) 'C change'],[char(0176) 'C change'],'change','% change','% change', ...
    '','','','','% change','% change','days change','days change'}
n.vars = length(vars);

% load AMM7 mesh
clear M
fname_mesh_h = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/ECOWINGS/RECICLE mesh files/mesh_hgr.nc';
fname_mask = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/ECOWINGS/RECICLE mesh files/mask.nc';
fname_mesh_z = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/ECOWINGS/RECICLE mesh files/mesh_zgr.nc';
% read T grid. note lon doesn't need wrapping.
M.lon_T = ncread(fname_mesh_h,'glamt');  % lon at T points
M.lat_T = ncread(fname_mesh_h,'gphit');  % lat at T point
n.lons = size(M.lon_T,1);
n.lats = size(M.lat_T,2);
fname_bathy = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/ECOWINGS/RECICLE mesh files/bathy_meter.nc';
M.bathymetry = ncread(fname_bathy,'Bathymetry');

% load averaging areas mask
fname = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/MATLAB/ECOWINGS/data/NWS_AMM7_box_mask15.nc';
M.boxmask = ncread(fname,'boxmask');
region_names = {'southern_North_Sea','CW_North_Sea','NW_North_Sea','English_Channel','Skaggerrak_Kattegat',...
    'Norwegian_Trench','Shetland_Shelf','Irish_Shelf','Irish_Sea','Celtic_Sea','Amorican_Shelf',...
    'offshelf_south','offshelf_north','CE_North_Sea','NE_North_Sea'};
M.xvec = M.lon_T(:,1);
M.yvec = M.lat_T(1,:)';

region_names = {'S_North_Sea','CW_North_Sea','NW_North_Sea','English_Channel','Skaggerrak_Kattegat',...
    'Norwegian_Trench','Shetland_Shelf','Irish_Shelf','Irish_Sea','Celtic_Sea','Amorican_Shelf',...
    'S_offshelf','N_offshelf','CE_North_Sea','NE_North_Sea'};

% define larger areas
aggregated.northern_North_Sea = [3 15];
aggregated.central_North_Sea = [2 14];
aggregated.North_Sea = [2 3 14 15];
aggregated.offshelf = [12 13];
aggregated.onshelf = setdiff([1:15],aggregated.offshelf);

larger_regions = fieldnames(aggregated);

all_regions = [region_names' ; larger_regions];
all_regions_long = {'southern North Sea','CW North Sea','NW North Sea','English Channel','Skaggerrak Kattegat',...
    'Norwegian Trench','Shetland Shelf','Irish Shelf','Irish Sea','Celtic Sea','Amorican Shelf',...
    'offshelf, south','offshelf, north','CE North_Sea','NE North_Sea', ...
    'northern North Sea','central North Sea', 'North Sea','offshelf','onshelf'};

n.regions = length(all_regions);

% time periods
time_periods = {'now','NF','FF'};
time_periods_long = {'Historic','Near future (2040s)','Far future (2090s)'};
n.tp=length(time_periods);

%% retrieve additive salinity Deltas
clear Delta_add

for m=1:n.mdls
    load(fullfile(infpath,['Delta_additive_SSS_' mdlName{m} '.mat']));
    Delta_add.(mdlName{m})= Delta_add_SSS;
end

%% get mean Delta_all over year
% note don't need to do this for bloom phenology

clear Delta_annual

for m=1:n.mdls
    for v=1:(n.vars-2)
        for t=2:n.tp
            % mean over months
            if strcmp(vars{v},'sss') % use additive for salinity
                Delta_annual.(mdlName{m}).(vars{v}).(time_periods{t}) = ...
                    mean(Delta_add.(mdlName{m}).(time_periods{t}),3,'omitmissing');
            else
                Delta_annual.(mdlName{m}).(vars{v}).(time_periods{t}) = ...
                    mean(Delta_all.(mdlName{m}).(vars{v}).(time_periods{t}),3,'omitmissing');
            end
        end
    end
end

for m=1:n.mdls
    for v=(n.vars-1):(n.vars)
        for t=2:n.tp
            if ~strcmp(mdlName{m},'HADGEM')
                % grab additive bloom Deltas
                Delta_annual.(mdlName{m}).(vars{v}).(time_periods{t}) = ...
                    Delta_bloom_add.(mdlName{m}).(vars{v}).(time_periods{t});
            end
        end
    end
end


%% get mean Delta by area

% do this for the original areas first
clear mean_Delta

for r=1:n.regions
    if r<=15
        idx=(M.boxmask==r);
    else
        idx=ismember(M.boxmask,aggregated.(all_regions{r}));
    end
    for v=1:n.vars
        for m=1:n.mdls
            if ~( strcmp(mdlName{m},'HADGEM') & (strcmp(vars{v},'bloom_start') | strcmp(vars{v},'bloom_duration') ))
                for t=2:n.tp
                    mean_Delta.(mdlName{m}).(vars{v}).(time_periods{t}).(all_regions{r}) = ...
                        mean(Delta_annual.(mdlName{m}).(vars{v}).(time_periods{t})(idx),'omitmissing');
                end
            end
        end
    end
end


%% as above but want SST and bloom phenology additive
% and the rest as % change

fileID = fopen('C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\RF zooplankton climate change paper\tables\predictors_projns_Delta__mod_SSS.csv','w');

fprintf(fileID,'%18s,%6s','Predictor variable','Region');
% write headers
fprintf(fileID,',%9s,%9s,%7s,%7s,%7s,%7s', ...
    'HADGEM_NF','IPSL_NF','GFDL_NF','HADGEM_FF','IPSL_FF','GFDL_FF');
fprintf(fileID,'\n');

for v=[1 3 4 5 10 11 12 13]%1:n.vars
    for r=[12 13 20] %1:length(all_regions)
        %disp(all_regions{r})
        fprintf(fileID,['%' num2str(length(vars{v})) 's'],vars{v});
        fprintf(fileID,[',%' num2str(length(all_regions{r})) 's'],all_regions{r});
        for t=2:n.tp
            for m=1:n.mdls

                % write SST & bloom dynamics as absolute change
                if strcmp(vars{v},'sst') | strcmp(vars{v},'sss') | (~strcmp(mdlName{m},'HADGEM') & (strcmp(vars{v},'bloom_start') | strcmp(vars{v},'bloom_duration')))
                    toWrite =  mean_Delta.(mdlName{m}).(vars{v}).(time_periods{t}).(all_regions{r});
                    fprintf(fileID,',%4.2f',toWrite);
                else
                    % write other variables as % change
                    if ~( strcmp(mdlName{m},'HADGEM') & (strcmp(vars{v},'bloom_start') | strcmp(vars{v},'bloom_duration') ))
                        toWrite = 100.0*(mean_Delta.(mdlName{m}).(vars{v}).(time_periods{t}).(all_regions{r})-1);
                        fprintf(fileID,',%4.1f',toWrite);
                    else
                        fprintf(fileID,',');
                    end
                end
            end
        end
        fprintf(fileID,'\n');
    end
end
fclose(fileID);

%% make plot
% easiest way is to load the table

delta_tbl = readtable('C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\RF zooplankton climate change paper\tables\predictors_projns_Delta__mod_SSS.csv');
% separate bar for each variable?
close all
figure
t=tiledlayout(8,2,'TileIndexing','columnmajor')

plot_cols = {3:5,6:8}
subtitles={'(a)','(b)'}
for i=1:2
    for v=1:n.vars
        idx = find(strcmp(delta_tbl.PredictorVariable,vars{v}));
        if ~isempty(idx)
            nexttile
            barh(delta_tbl{idx,plot_cols{i}},'LineStyle','none')
            if v==1 
                title({subtitles{i} , vars_long{v}})
            else
            title(vars_long{v})
            end
            xlabel(vars_delta_units{v})
            if i==1
                yticklabels(all_regions_long([12 13 20]))
            else
                yticklabels([])
            end
            if v==1 & i==1
                lg=legend(mdlName_long,'Orientation','horizontal');
            end

            xlim_min = min(delta_tbl{idx,3:8},[],'all');
            xlim_max = max(delta_tbl{idx,3:8},[],'all');

            xlim([xlim_min xlim_max])

        end
    end
end
lg.Layout.Tile = 'south';
t.Padding='compact';
t.TileSpacing="compact";

plotname='climate_projns_deltas.png'
print(fullfile(plotpath,plotname),'-dpng')




