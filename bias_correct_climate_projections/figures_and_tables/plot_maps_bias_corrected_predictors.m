%% make contour plots for decadal annual predictors
clear, close all, clc
mdlName = {'HADGEM','RECICLE_IPSL','RECICLE_GFDL'};
mdlName_long = {'HADGEM','IPSL','GFDL'};
n.mdls=length(mdlName);

plotpath = 'C:\Users\xqb21120/OneDrive - University of Strathclyde\Documents\RF zooplankton climate change paper\plots';

infpath ='C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\ECOWINGS\climate projections bias correction';

% load predictors by model
% e.g. sst lon*lat*month*year
clear predictors_all
for m=1:n.mdls
    disp(['loading ' mdlName{m} ' predictors'])
    load(fullfile(infpath,['bias_corrected_' mdlName{m} '_predictors_incl_no3.mat'])); % should save this on T7 - is too large for harddrive
    predictors_all.(mdlName{m})=predictors;
end

vars = fieldnames(predictors);
n.vars = length(vars);
vars_long = {'sea surface temperature','salinity', ...
    'mixed layer depth','monthly NO_3','chlorophyll','winter nitrate','mo SST anom','mo NO_3 anom','bloom start','bloom duration'};

 
time_periods = {'now','NF','FF'};
tp_long={'recent past','near future','far future'};
time_decades = [2010,2040,2090];
n.tp = length(time_periods);

%% load AMM7 mesh - can use same for all projections because only using lat
% lon, not depth
clear M

% same for both RECICLE runs
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

%%
close all
    figure
    t=tiledlayout(3,2,'TileIndexing','rowmajor');
for v=[1:2 6]%[1:8 10:n.vars] % don't plot SPG proxy


    for i=[1 3]
    for m=2%:n.mdls

        if strcmp(mdlName{m},'HADGEM') & (strcmp(vars{v},'bloom_start') | strcmp(vars{v},'bloom_duration'))
            disp('no bloom dynamics')
        else
            if strcmp(vars{v},'bloom_start') | strcmp(vars{v},'bloom_duration')
                data = ( mean( predictors_all.(mdlName{m}).(vars{v}).(time_periods{i}), 3) );
            else
                data = ( mean( predictors_all.(mdlName{m}).(vars{v}).(time_periods{i}), [3 4]) );
            end

            if i==1 & m>2
                disp('not plotting')
            else
            nexttile

            if v==1
                contourf(M.lon_T,M.lat_T,data,4:2:22,'ShowText','off','LineStyle','none');
                caxis([4 22])
                cmocean('thermal',10)
            cb=colorbar;
            title(cb,'deg C')
            elseif v==2
                contourf(M.lon_T,M.lat_T,data,32:0.5:36,'ShowText','off','LineStyle','none');
                caxis([32 36])
                cmocean('haline',8)
                            cb=colorbar;

            else
                [c,h]=contourf(M.lon_T,M.lat_T,data,0:2:16,'ShowText','off','LineStyle','none') ;
                caxis([0 16])
                cmocean('matter',8)
                            cb=colorbar;
                            title(cb,'mmol.m^{-3}')
                           
clabel(c,h,'Color','w')

            end
            if i==1 & v==1
                title('Recent past')            
                %title(['Recent past ' vars_long{v}],'Interpreter','none')

            elseif v==1
                title('2090-2100')
            %title([tp_long{i} ' ' mdlName_long{m} ' ' vars_long{v}],'Interpreter','none')
            end
            if i==1
                ylabel(vars_long{v})
            end
            xticklabels([])
            yticklabels([])
            end
        end
    end
    end
end
t.Padding="compact";
t.TileSpacing="compact";

%print(fullfile(plotpath,'projected_sst_sss_NO3'),'-dpng')
print(fullfile('C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\ECOWINGS\topic sheets','projected_sst_sss_NO3'),'-dpng')
