%% high level file for processing the random forest predictions
% load each time step and collate: function process_RF_predictions.m

% amended Oct 2024 for R output with all Z classes in one file
% amended Feb 2025 to include boxmask
% amended Apr 2025 to generalise to any Z classes
% modified Dec 2025 for manuscript

clear, close all, clc
fpath='C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/MATLAB/ECOWINGS/randomForest/output';

plotpath = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/MATLAB/ECOWINGS/plots';

% which projection?
oceanModelName={"HADGEM","RECICLE_IPSL","RECICLE_GFDL"};
n.mdls = length(oceanModelName);

% which area RF trained on?
trainedOn = "NEAtlantic"; 
% which area predicted on?
predictedOn = "NEAtlantic"; 

 % time periods
decades = [2010 2040 2090];
n.dec=3;
dt="monthly";

%% loop over time periods
for m=2%:n.mdls

    if strcmp(oceanModelName{m},"HADGEM")
        mdlStr="reduced_with_no3";% "reduced_exclude_coastal"; 
    else
        mdlStr= "full_with_no3"; %"full_exclude_coastal";
    end
    disp("processing RF predictions for: "+mdlStr)

    disp("trained on: "+trainedOn)
    disp("predicted on: "+predictedOn)
    disp("for projn "+oceanModelName{m})

    for d = 1%:n.dec

        disp(['processing for ' num2str(decades(d)) 's'])
        year_start = decades(d);
        year_end = decades(d)+9;

        % collate time step files & save to file
        outname = "predicted_"+mdlStr+"_projn_"+oceanModelName{m}+"_"+trainedOn+"_"+dt+"_"+num2str(year_start)+"_"+num2str(year_end)+".mat"
        Z = process_rf_predictions_monthly(mdlStr,trainedOn,predictedOn,year_start,year_end,oceanModelName{m});

        % remove rows with NaT
        Z(isnat(Z.t),:)=[];

        % save to file
        save(fullfile(fpath,outname),"Z",'-v7.3');
    end
end

