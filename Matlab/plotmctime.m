function [h, mctime, mcaverages] = plotmctime(energies)
%PLOTMCTIME(ENERGIES)
%  ENERGIES:  vector of the energy samples


    % Choose a time step such that we have at most 1000 points.
    % Decomment the commented line below if you want to manually override.
    mctimestep = length(energies)/1000;
    %mctimestep = 100;

    % Do a preallocation of the vector for speed
    iteration = 1;
    mcaverages = zeros(1,length(energies)/mctimestep);

    % Generate the MC time points
    mctime = mctimestep:mctimestep:length(energies);

    % Do the calcs
    for timeint = mctime
        mcaverages(iteration) = mean(energies(1:timeint));
        iteration = iteration + 1;
    end

    % Create figure
    MCTimeFigure = figure('PaperOrientation','landscape','PaperType','A3');

    % Create axes
    axes1 = axes('Parent',MCTimeFigure);
    hold(axes1,'on');

    % Create plot
    h = plot(mctime, mcaverages, 'LineWidth',2);

    % Create xlabel
    xlabel('Monte Carlo steps','FontSize',24);

    % Create title
    title('$\langle V \rangle$ vs. Number of Monte Carlo steps','FontSize',24,...
        'Interpreter','latex');

    % Create ylabel
    ylabel('$\langle V \rangle$ [$\varepsilon$]','FontSize',24,...
        'Interpreter','latex');

    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'FontSize',18);

end