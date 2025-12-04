%% plot projection monthly Deltas

clear, close all, clc
mdlName = {'HADGEM','RECICLE_IPSL','RECICLE_GFDL'};
mdlName_long = {'HADGEM','RECICLE IPSL','RECICLE GFDL'};
n.mdls=length(mdlName);

plotpath = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/RF zooplankton climate change paper/plots/';
infpath = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/MATLAB/ECOWINGS/climate projections bias correction/output';

% load Delta values
clear Delta_all
load(fullfile(infpath,'Delta_all_models_updated.mat'))

% load bloom phenology Deltas
for m=1:n.mdls
    if ~(strcmp(mdlName{m},'HADGEM'))
        load(fullfile(infpath,['Delta_bloom_updated_' mdlName{m} '.mat']));
        Delta_bloom.(mdlName{m})= Delta;
        Delta_all.(mdlName{m}).bloom_start = Delta.bloom_start;
        Delta_all.(mdlName{m}).bloom_duration = Delta.bloom_duration;
    end
end

clear Delta Delta_bloom

vars = fieldnames(Delta_all.RECICLE_GFDL);
n.vars = length(vars);

vars_long = {'sea surface temperature','near bottom temperature','sea surface salinity', ...
    'mixed layer depth','chlorophyll','','','','','','winter nitrate','bloom start','bloom duration'};

vars_units = {[char(0176) 'C'],[char(0176) 'C'],'','log_{10} m','mg.m^{-3}', ...
    '','','','','mmol.m^{-3}','mmol.m^{-3}','yearday','days'};

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

% time periods to plot
time_periods = {'now','NF','FF'};
time_periods_long = {'Historic','Near future (2040s)','Far future (2090s)'}
n.tp=length(time_periods);

% plot Deltas
mo.annual = 1:12;
mo.summer = 3:8;
mo.winter = [1:2 9:12];

season = fields(mo);

%% make histograms to help with plot limits

close all
plotVars = [1 3:5 11:n.vars];
n.plotV=length(plotVars)
    figure
    t=tiledlayout(n.plotV,n.mdls);
    for p=[1 3:5 11:n.vars]
        for m=1:n.mdls
            nexttile
            for i=2:n.tp % plot NF and FF

            % plot annual mean unless is already an annual metric
            if ~(ismember(vars{p},{'bloom_start','bloom_duration'}) & strcmp(mdlName{m},'HADGEM'))
                if ismember(vars{p},{'no3_winter','bloom_start','bloom_duration'})
                    data = Delta_all.(mdlName{m}).(vars{p}).(time_periods{i})  ;
                else
                    data = mean(Delta_all.(mdlName{m}).(vars{p}).(time_periods{i}),3);
                end
                  histogram(data,'LineStyle','none');
                  data_mean = mean(data,'all','omitmissing');
                  data_sd = std(data(:),'omitmissing');

                  if i==1
                    xlim([data_mean-2*data_sd,data_mean+2*data_sd])
                  else
                      xlim_old = xlim;

                      xlim([min(data_mean-3*data_sd,xlim_old(1)), ...
                          max(data_mean+3*data_sd,xlim_old(2))])
                  end
                  
                hold on
            end
            end
            if p==1
                title(mdlName_long{m})
            end
            
            %if m==1
                xlabel([vars_long{p} ' (' vars_units{p} ')'],'fontsize',8,'rotation',0)
            %end
            if strcmp(vars{p},'sst')
                      xline(0)
                  else
                      xline(1)
                  end
        end
    end
    t.Padding="tight";
    t.TileSpacing="tight";
    lg=legend('Near future (2040s)','Far future (2090s)','Orientation','horizontal');
    lg.Layout.Tile = 'North';
plotname = ['hist_Delta_predictors'];
sgtitle('Projected change in each variable (additive for SST, multiplicative otherwise)')
print(fullfile(plotpath,plotname),'-dpng')

%% plot Deltas for paper Supp Info
close all

plotLims.sst = [-1 4];
plotLims.sss = [0.9 1.05];
plotLims.mld = [0.5 1.5];
plotLims.chl = [0.5 3];
plotLims.no3_winter = [0.2 1.5];
plotLims.bloom_start = [0.8 1.2];
plotLims.bloom_duration = [0.8 1.3];

plotVars = [1 3:5 11:n.vars];
n.plotV=length(plotVars)
for i=2:n.tp % plot NF and FF
    figure
    t=tiledlayout(n.plotV,n.mdls);
    for p=[1 3:5 11:n.vars]

        for m=1:n.mdls

            nexttile
            % plot annual mean unless is already an annual metric
            if ~(ismember(vars{p},{'bloom_start','bloom_duration'}) & strcmp(mdlName{m},'HADGEM'))
                if ismember(vars{p},{'no3_winter','bloom_start','bloom_duration'})
                    pcolor(M.lon_T, M.lat_T,Delta_all.(mdlName{m}).(vars{p}).(time_periods{i})), shading flat
                else
                    pcolor(M.lon_T, M.lat_T,(mean(Delta_all.(mdlName{m}).(vars{p}).(time_periods{i}),3))), shading flat
                end
            end

            if p==1
                title(mdlName_long{m})
            end

            
            if m==1
                ylabel({vars_long{p},vars_units{p}},'fontsize',6)
            end

            caxis(plotLims.(vars{p}))

            cb=colorbar;
            %title(cb,['/Delta' vars{p}]);
            if strcmp(vars{p},'sst')
                cmocean('balance',10,'pivot',0)
            else
                cmocean('balance',10,'pivot',1)
            end

        end
    end
    t.TileSpacing='tight';
t.Padding='tight';
plotname = ['Delta_predictors_',time_periods{i}];
print(fullfile(plotpath,plotname),'-dpng')
end

%% load SPG proxy projections
fdir = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\ECOWINGS\climate projections bias correction\'
clear SPG_proxy
for m=1:n.mdls
    fname = fullfile(fdir,['bias_corrected_' mdlName{m} '_predictors.mat']);
    disp(['loading predictors: ' mdlName{m}])
    load(fname)
    SPG_proxy.(mdlName{m}) = predictors.SPG_proxy_anom
end

%% plot SPG proxy projections
close all
figure
plot(1:14,SPG_proxy.(mdlName{1}).now,'-k','linewidth',1.5,'DisplayName','Historical (2010-2020)')
hold on

lines_type = {'--',':'}
lines_col = {'r','g','b'}
for m=1:n.mdls
    for i=2:n.tp
        plot(1:14,SPG_proxy.(mdlName{m}).(time_periods{i}), ...
            'DisplayName',[mdlName{m} ' ' time_periods_long{i}], ...
            'LineStyle',lines_type{i-1},'Color',lines_col{m})
    end
    grid on


end 

xlabel('year')
ylabel('upper ocean salinity anomaly')
legend('location','eastoutside','Interpreter','none')
ylim([-0.35 0.1])
xlim([1 10])

print(fullfile(plotpath,'SPG_projns_all_models.png'),'-dpng')

