function plotautocorrelation(data, varargin)
%PLOTAUTOCORRELATION(DATA, CUTOFF)
%  DATA:  Vector of data
%  CUTOFF:  Cutoff of the data in the calculation of the autocorrelation.
%           Optional argument, default is CUTOFF=100000.

    
    % Parse input
    iP = inputParser;
    iP.FunctionName = 'PLOTAUTOCORRELATION';
    iP.addRequired('data',@isnumeric);
    iP.addOptional('cutoff',100000,@isnumeric);    
    iP.parse(data, varargin{:});
    cutoff = iP.Results.cutoff;

    % Create figure
    AutocorrelationFigure = figure('PaperOrientation','landscape','PaperType','A3');

    % Create axes
    AutocorrelationAxes = axes('Parent',AutocorrelationFigure);
    hold(AutocorrelationAxes,'on');

    % Create plot
    plot(autocorr(data,cutoff),'LineWidth',2);

    % Create xlabel
    xlabel('$\tau$','FontSize',24,'Interpreter','latex');

    % Create title
    title('Autocorrelation of $V(t)$','FontSize',24,'Interpreter','latex');

    % Create ylabel
    ylabel('$\frac{\Big\langle \big(V(t) - \langle V \rangle \big)\big(V(t+\tau) - \langle V \rangle \big) \Big\rangle}{\Big\langle \big( V(t) - \langle V \rangle \big)^2 \Big\rangle}$',...
        'FontSize',24,...
        'Interpreter','latex');

    % Set the X-limits of the axes
    xlim(AutocorrelationAxes,[0 cutoff]);
    
    % Activate the box
    box(AutocorrelationAxes,'on');
    
    % Set the remaining axes properties
    set(AutocorrelationAxes,'FontSize',18);

end