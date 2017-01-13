function  [ Vmean, Verror, Cv, Cverror, Pmean, Perror ] = IsothermAnalysis( NumPart, dataLines, TemperatureRange, DensityRange )
%MULTIPLETEMPERATUREANALYSIS Analyze multiple datasets produced by LJMC
% WARNING: Since I'm lazy, NumPart and Dimension are numbers.
% TemperatureRange and MCStepRange are vectors. The other vars are STRINGS.
%
%
% Example:
%   Write an example;


    %% Do preallocations
    Vmean = zeros(length(TemperatureRange),length(DensityRange));
    Verror = zeros(length(TemperatureRange),length(DensityRange));
    Cv = zeros(length(TemperatureRange),length(DensityRange));
    Cverror = zeros(length(TemperatureRange),length(DensityRange));
    Pmean = zeros(length(TemperatureRange),length(DensityRange));
    Perror = zeros(length(TemperatureRange),length(DensityRange));
    totalSets = length(TemperatureRange)*length(DensityRange);

    
    %% Import the data and do the things
    for indexT = 1:length(TemperatureRange)
        for indexRHO = 1:length(DensityRange)
            % Import data in a faster way
            confString = strcat(...
                            '[T',num2str(TemperatureRange(indexT),'%.2f'),'-',...
                            'RHO',num2str(DensityRange(indexRHO),'%.2f'),...
                            ']');
%             itemString = strcat(...
%                             'T',num2str(TemperatureRange(indexT)*100,'%d'),...
%                             'RHO',num2str(DensityRange(indexRHO)*100,'%d'));
                        
            data = ImportIsotherm(strcat('out_file_-_',confString), dataLines);
                        
            % Calculate Vmean and its error
            Vmean(indexT,indexRHO) = mean(data(:,1)/NumPart);
            Verror(indexT,indexRHO) = binanalysis(data(:,1)/NumPart, 'No');
            
            % Calculate Cv and its error
            [Cv(indexT,indexRHO), Cverror(indexT,indexRHO)] = specificheat(data(:,1), TemperatureRange(indexT), NumPart);

            % Calvulate Pmean and its error
            Pmean(indexT,indexRHO) = mean(data(:,2));
            Perror(indexT,indexRHO) = binanalysis(data(:,2), 'No');
            
            % Save only the data and clear it afterwards
%             save(strcat('dataWorkspace_-_',confString,'.mat'), '-struct', 'data', itemString);
%             data.(itemString) = [];
            save(strcat('dataWorkspace_-_',confString,'.mat'), 'data');
            clear data;
            
            % Backup the current workspace
            save('workspace.mat',...
                'Vmean','Verror',...
                'Cv','Cverror',...
                'Pmean','Perror');
            
            % Print out some infos
            fprintf( 'Done dataset %d of %d.\n', (indexT-1)*length(DensityRange)+indexRHO, totalSets);
            
        end
    end
    

end

