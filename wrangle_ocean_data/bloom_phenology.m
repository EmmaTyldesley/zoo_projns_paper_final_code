% Get phytoplankton bloom phenology using daily AMM7 reanalysis P
% concentration. Bloom start and end defined as days with 15% and 85%
% cumulative annual total respectively. Duration is number of days between
% start and end. Metrics calculated annually by lat lon grid-cell.

% EJT tidied Dec 2025 for manuscript submission

clear, close all

saveData=0;

% get bathymetry - from Copernicus reanalysis - see README.txt
filepath = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\AMM7 masks\';
filename_bathy = [filepath 'NWS-MFC_004_001_mask_bathy.nc'];
M.deptho = ncread(filename_bathy,'deptho');
M.lon = ncread(filename_bathy,'longitude');
M.lat = ncread(filename_bathy,'latitude');

% phytoplankton concentration from Copernicus reanalysis - see README.txt
fpath = 'D:\NWS_multiyear\bgc\phyc\daily\';

yrs = 1998:2022;
n.yrs = length(yrs);

blank = nan(297,375,n.yrs); % lon x lat x year
bloom.start = blank;
bloom.duration = blank;
bloom.stop = blank;
bloom.years = yrs;

%% calculate bloom dynamics
% start = day of year when cumulative total is 15% annual total
% duration = number of days between 15% and 85% annual total
for i = 1:n.yrs
    yr=yrs(i);
    Pconc_annual = [];
    t_annual = [];
    yrStr = num2str(yr)

    % get filenames
    for mo=1:12
        if mo<10
            moStr = ['0' num2str(mo)]
        else
            moStr = num2str(mo)
        end

        flist = dir([fpath yrStr '\' moStr]);

        for f=3:length(flist)
            fname = fullfile(flist(f).folder,flist(f).name)

            % read variables
            if ~isfield(M,'depth')
                disp('getting grid info...')
                M.depth = ncread(fname,"depth");
                disp('done.')
            end

            M.Pconc = ncread(fname,"phyc");
            M.t = ncread(fname,"time");
            M.t = datetime(M.t,'ConvertFrom','epochtime','Epoch','01-Jan-1970 00:00:00');

            % integrate P over top 50 m => P 297x375
            zlayers=find(M.depth<=50); %
            M.Pconc_int = integrateByDepth(M.depth(zlayers),M.depth(zlayers(1)),M.depth(zlayers(end)),M.Pconc); %integrateByDepth(zlayers,ztop,zbot,var)

            % stack P(t) to make P_all 297x375x365(6)
            Pconc_annual = cat(3,Pconc_annual,M.Pconc_int);
            t_annual = [t_annual; M.t];
        end
    end

    % get bloom metrics for one year
    Pcumsum = cumsum(Pconc_annual,3);      % get cumulative sum along t dimension
    start_conc = 0.15*Pcumsum(:,:,end);    % P conc for start of bloom
    end_conc = 0.85*Pcumsum(:,:,end);      % P conc for end of bloom

    for j=1:length(M.lon)
        for k=1:length(M.lat)
            idx = find( Pcumsum(j,k,:)>=start_conc(j,k) & Pcumsum(j,k,:)<=end_conc(j,k) );
            if ~isempty(idx)
                bloom.start(j,k,i) = day(t_annual(idx(1)),'dayofyear');
                bloom.stop(j,k,i) = day(t_annual(idx(end)),'dayofyear');
                bloom.duration(j,k,i) = length(idx); % days
            end
        end
    end
end

% save
if saveData
    save(['data/bloom_1998_2022.mat'],'bloom')
end

