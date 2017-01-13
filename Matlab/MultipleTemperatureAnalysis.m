function  [ data, Vmean, Verror, Cv, Cverror ] = MultipleTemperatureAnalysis( NumPart, Dimension, TemperatureRange, MCNumSteps, MCStepRange, Thermalization )
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

    
    %% Import the data
    parfor index = 1:length(TemperatureRange)
        data{index} = importfile(strcat('out_file_-_[N',num2str(NumPart),...
            '_-_D',num2str(Dimension),'_-_T',num2str(TemperatureRange(index),'%.1f'),...
            '_-_MCS',MCNumSteps,'_-_S',num2str(MCStepRange(index),'%.2f'),...
            '_-_THS',Thermalization,']'));
    end
    
    
    %% Do all the calcs
    % Calvulate Vmean and its error
    parfor index = 1:length(TemperatureRange)
        Vmean(index) = mean(data{index});
        Verror(index) = binanalysis(data{index});
    end
    
    % Calculate Cv and its error
    parfor index = 1:length(TemperatureRange)
        [Cv(index), Cverror(index)] = specificheat(data{index}, TemperatureRange(index));
    end
    
    
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

