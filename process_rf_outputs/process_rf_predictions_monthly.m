function [Z] = process_rf_predictions_monthly(mdlStr,trainedOn,predictedOn,year_start,year_end,oceanModelName)
    % process RF predictions by zooplankton class
    % load timesteps and collate into struct
    % rf predictions carried out in R
    

    % - amended Oct 2024 for all Z classes in one R output file
    % - amended Apr 2025 to read Z classes from example file (move beyond
    % fixed taxa)
    % - renamed Dec 2025 for manuscript
    
    infpath = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/temp_rf_output_files'; % rf output on amm7 reanalysis
    outfpath = 'C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/MATLAB/ECOWINGS/DATA';
        
    yrs = year_start:year_end;
    n.yrs = length(yrs);
    
    % --- Make N by 2 matrix of fieldname + data type
    % read one file to check size
    % example filename: predicted_ZE_full_with_no3_NEAtlantic_NEAtlantic_RECICLE_GFDL_monthly_201001

    % note these files are stored temporarily in folder
    % temp_rf_output_files and removed after processing
    fname = fullfile(infpath,'predicted_ZE_'+mdlStr+'_'+trainedOn+'_'+predictedOn+'_'+oceanModelName+'_monthly_'+num2str(yrs(1))+'01.txt');
    disp(['predicted_ZE_'+mdlStr+'_'+trainedOn+'_'+predictedOn+'_'+oceanModelName+'_monthly_'+num2str(yrs(1))+'01.txt'])
    test = readtable(fname);
    n.rows = height(test);
    taxa = convertCharsToStrings(test.Properties.VariableNames(6:end))';
    n.taxa = height(taxa);
    variables_taxa = [taxa repmat("double",n.taxa,1)];
    variable_names_types = [ ["t", "datetime"]; ...
        ["x", "double"]; ...
        ["y", "double"]; ...
        ["depth", "double"];...
        ["boxmask", "double"]; ...
        variables_taxa ];
 
    % initialise output struct
    predicted_zed_all = table('Size',[n.rows*366*n.yrs height(variable_names_types)],...
        'VariableNames', variable_names_types(:,1),...
        'VariableTypes', variable_names_types(:,2));

    % loop over timesteps
    rownum=0;   % reset row number
    for yr = yrs
        YYYY = num2str(yr);
        for mo=1:12
            MM=sprintf('%02d',mo);

            fname = fullfile(infpath,'predicted_ZE_'+mdlStr+'_'+trainedOn+'_'+predictedOn+'_'+oceanModelName+'_monthly_'+YYYY+MM+'.txt');
            disp(['predicted_ZE_'+mdlStr+'_'+trainedOn+'_'+predictedOn+'_'+oceanModelName+'_monthly_'+YYYY+MM+'.txt'])
            if isfile(fname) & ~isempty(readtable(fname))
                predicted_zed = readtable(fname);
                predicted_zed_all( (rownum+1):(rownum+height(predicted_zed)) , :) = predicted_zed;
                rownum = rownum+height(predicted_zed);
            else
                error(['warning! no file or empty file for ' YYYY MM])
            end
        end
    end
    Z = predicted_zed_all; 
end


