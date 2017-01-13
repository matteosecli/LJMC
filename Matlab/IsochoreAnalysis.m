function  [ data, Vmean, Verror, Cv, Cverror, Pmean, Perror ] = IsochoreAnalysis( NumPart, TemperatureRange )
%MULTIPLETEMPERATUREANALYSIS Analyze multiple datasets produced by LJMC
% WARNING: Since I'm lazy, NumPart and Dimension are numbers.
% TemperatureRange and MCStepRange are vectors. The other vars are STRINGS.
%
%
% Example:
%   Write an example;

    %% Do preallocations
    data = cell(1,length(TemperatureRange));
    Vmean = zeros(1,length(TemperatureRange));
    Verror = zeros(1,length(TemperatureRange));
    Cv = zeros(1,length(TemperatureRange));
    Cverror = zeros(1,length(TemperatureRange));
    Pmean = zeros(1,length(TemperatureRange));
    Perror = zeros(1,length(TemperatureRange));

    
    %% Import the data
    parfor index = 1:length(TemperatureRange)
        data{index} = importdata(strcat('out_file_-_[T',...
                        num2str(TemperatureRange(index),'%.1f'),']'));
    end
    
    save workspace
    
    %% Do all the calcs
    % Calvulate Vmean and its error
    parfor index = 1:length(TemperatureRange)
        Vmean(index) = mean(data{index}(:,1)/NumPart);
        Verror(index) = binanalysis(data{index}(:,1)/NumPart, 'No');
    end
    
    % Calculate Cv and its error
    parfor index = 1:length(TemperatureRange)
        [Cv(index), Cverror(index)] = specificheat(data{index}(:,1), TemperatureRange(index), NumPart);
    end
    
    % Calvulate Pmean and its error
    parfor index = 1:length(TemperatureRange)
        Pmean(index) = mean(data{index}(:,2));
        Perror(index) = binanalysis(data{index}(:,2), 'No');
    end
    
    save workspace
    
    
    %% Plot the temperature
    % Create figure
    figure1 = figure;
    
    % Create axes
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    
    % Create errorbar
    errorbar(TemperatureRange,Cv,Cverror,'MarkerSize',8,'LineStyle','none','LineWidth',2);
    
    % Create xlabel
    xlabel('T [\epsilon/k_B]','FontSize',22);
    
    % Create title
    title('C_V vs T','FontWeight','bold','FontSize',22);
    
    % Create ylabel
    ylabel('C_V [k_B]','FontSize',22);
    
    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'FontSize',18); 
    
    % Print out some infos
    %fprintf( '\n V = %.4f +/- %.4f \n Cv = %.4f +/- %.4f \n\n', Vmean, Verror, Cv, Cverror);

end

