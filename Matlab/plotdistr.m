function [h meanvalue stdev] = plotdistr(energies)
%PLOTDISTR Plot the energy distribution
%   Plot the distribution of the LJ potential values, and returns
%   [h mean std] where 'h' is the handle of the histogram, 'mean' is the
%   mean value of the potential, and 'std' is its standard deviation


    % Create figure
    HistogramFigure = figure('PaperOrientation','landscape','PaperType','A3');

    % Create axes
    HistogramAxes = axes('Parent',HistogramFigure);

    % Create histogram
    h = histogram(energies,'DisplayName','data','Parent',HistogramAxes,...
        'EdgeColor',[0 0.447058823529412 0.741176470588235],...
        'NumBins',1000);

    % Create xlabel
    xlabel('V [$\varepsilon$]','FontSize',24,'Interpreter','latex');

    % Create ylabel
    ylabel('Countings','FontSize',24,'Interpreter','latex');

    box(HistogramAxes,'on');
    % Set the remaining axes properties
    set(HistogramAxes,'FontSize',18);
    
    % Calculate some useful statistical values
    meanvalue = mean(energies);
    stdev = std(energies);
    
end