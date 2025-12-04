%% combine Delta values calculated separately for phys, chl and bloom metrics

infpath = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\ECOWINGS\climate projections bias correction';
outfpath = 'C:\Users\xqb21120\OneDrive - University of Strathclyde\Documents\MATLAB\ECOWINGS\climate projections bias correction\output';

mdlName = {'HADGEM','RECICLE_IPSL','RECICLE_GFDL'}
n.mdls = length(mdlName);

% collate Delta values
clear Delta_all
for m=1:n.mdls
    disp(mdlName{m})
    % load Deltas
    load(fullfile(infpath,['Delta_SST_NBT_SSS_MLD_updated_' mdlName{m} '.mat'])); % Delta.<var>.<NF/FF> lon*lat*12
    Delta_all.(mdlName{m}) = Delta;
    load(fullfile(infpath,['Delta_CHL_PFT_NO3_' mdlName{m} '.mat']));
    P=fieldnames(Delta);
    for l=1:length(P)
        Delta_all.(mdlName{m}).(P{l}) = Delta.(P{l});
    end
    clear Delta
end

% save to file
fname = fullfile(outfpath,'Delta_all_models_updated.mat');
disp(fname)
save(fname,"Delta_all")
