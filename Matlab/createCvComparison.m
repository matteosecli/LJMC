function createCvComparison(TemperatureRange, Vmean, Verror, Cv, Cverror, Name)
%CREATEFIT(TEMPERATURERANGE,VMEAN3,VWEIGHTS3)
%  Create a fit.
%
%  Data for 'FitVrho3' fit:
%      X Input : TemperatureRange
%      Y Output: Vmean3
%      Weights : Vweights3
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.


% Define weights
Vweights = 1./Verror.^2;

%% Fit: 'FitVrho3'.
[xData, yData, weights] = prepareCurveData( TemperatureRange, Vmean, Vweights );

% Set up fittype and options.
ft = fittype( 'poly4' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Weights = weights;

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'FitVrho3' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Vmean3 vs. TemperatureRange with Vweights3', 'FitVrho3', 'Location', 'NorthEast' );
% % Label axes
% xlabel TemperatureRange
% ylabel Vmean3
% grid on

% Calculate the derivative
xDummy = linspace(min(xData), max(xData), 1000);
Dfitresult = differentiate(fitresult, xDummy);


%% Do the plot

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create the Dfit line
line1 = plot(xDummy, Dfitresult,...
    'Parent',axes1,...
    'DisplayName','Derivative of the energy fit, 3rd degree poly',...
    'Marker','none',...
    'LineStyle','-',...
    'LineWidth',2,...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625]);

% Create the Cv errorbar
errorbar1 = errorbar(TemperatureRange,Cv,Cverror,...
    'Parent',axes1,...
    'DisplayName',Name,...
    'MarkerSize',8,...
    'LineStyle','none',...
    'LineWidth',2,...
    'Color',[0 0.447058826684952 0.74117648601532]);

% Create xlabel
xlabel('$T$ [$\varepsilon/k_B$]','FontSize',22,'Interpreter','latex');

% Create title
title('$C_V$ vs $T$','FontWeight','bold','FontSize',24,...
    'Interpreter','latex');

% Create ylabel
ylabel('$C_V$ [$\varepsilon k_B^{-3}$]','FontSize',22,'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',18);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','best','Interpreter','latex','FontSize',18);

    

