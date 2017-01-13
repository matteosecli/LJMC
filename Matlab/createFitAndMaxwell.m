function [fitresult, fitZeros, TvsRho, TvsRhoErrors, MaxValue, MaxValueError] = createFitAndMaxwell(TemperatureRange, DensityRange, Pmean, Perror, Psat, JohnsonParams, JohnsonData)
%CREATEFIT(DENSITYRANGE,TESTY,TESTW)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : DensityRange
%      Y Output: testY
%      Weights : testW
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.



%% Initialization

% Variables
fitresult = cell(1,length(TemperatureRange));
gof = cell(1,length(TemperatureRange));
legendVector = zeros(1,length(TemperatureRange));
syms f(x);
syms f(rho);
fitZeros = cell(1,length(TemperatureRange));

% Do the fucking crazy fit
crazyfitresult = createCrazyFit(TemperatureRange, DensityRange, Pmean-Psat');

% Create figure
figure1 = figure('PaperOrientation','landscape','PaperType','A3');

% Create axes
axes1 = axes('Parent', figure1, 'FontSize', 18);
hold(axes1,'on');


%% Fit for each temperature
for indexT = 1:length(TemperatureRange)
    Temperature = TemperatureRange(indexT);
    [xData, yData, weights] = prepareCurveData( DensityRange, Pmean(indexT,:), 1./Perror(indexT,:).^2 );
    
    T = Temperature;    % It's stupid, but this way I don't have to change the whole fucking function
    f(rho) = rho.*T+rho.^2.*(crazyfitresult.x1.*T+crazyfitresult.x2.*sqrt(T)+crazyfitresult.x3+crazyfitresult.x4./T+crazyfitresult.x5./T.^2) ...
             + rho.^3.*(crazyfitresult.x6.*T+crazyfitresult.x7+crazyfitresult.x8./T+crazyfitresult.x9./T.^2) + rho.^4.*(crazyfitresult.x10.*T+crazyfitresult.x11+crazyfitresult.x12./T) ...
             + rho.^5.*(crazyfitresult.x13) + rho.^6.*(crazyfitresult.x14./T+crazyfitresult.x15./T.^2) + rho.^7.*(crazyfitresult.x16./T) ...
             + rho.^8.*(crazyfitresult.x17./T+crazyfitresult.x18./T.^2) + rho.^9.*(crazyfitresult.x19./T.^2) ...
             + ( rho.^3.*(crazyfitresult.x20./T.^2+crazyfitresult.x21./T.^3) ...
               + rho.^5.*(crazyfitresult.x22./T.^2+crazyfitresult.x23./T.^4)   + rho.^7.*(crazyfitresult.x24./T.^2+crazyfitresult.x25./T.^3) ...
               + rho.^9.*(crazyfitresult.x26./T.^2+crazyfitresult.x27./T.^4)   + rho.^11.*(crazyfitresult.x28./T.^2+crazyfitresult.x29./T.^3) ...
               + rho.^13.*(crazyfitresult.x30./T.^2+crazyfitresult.x31./T.^3+crazyfitresult.x32./T.^4) ) .* exp(-crazyfitresult.g.*rho.^2);

    % Set up fittype and options.
    ft = fittype( 'poly5' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf 0];
    opts.Upper = [Inf Inf Inf Inf Inf 0];
    opts.Weights = weights;

%     % Set up fittype and options.
%     ft = fittype( strcat('rho*',num2str(Temperature),'+rho^2*(x1*',num2str(Temperature),'+x2*sqrt(',num2str(Temperature),')+x3+x4/',num2str(Temperature),'+x5/',num2str(Temperature),'^2) + rho^3*(x6*',num2str(Temperature),'+x7+x8/',num2str(Temperature),'+x9/',num2str(Temperature),'^2) + rho^4*(x10*',num2str(Temperature),'+x11+x12/',num2str(Temperature),') + rho^5*(x13) + rho^6*(x14/',num2str(Temperature),'+x15/',num2str(Temperature),'^2) + rho^7*(x16/',num2str(Temperature),') + rho^8*(x17/',num2str(Temperature),'+x18/',num2str(Temperature),'^2) + rho^9*(x19/',num2str(Temperature),'^2) + ( rho^3*(x20/',num2str(Temperature),'^2+x21/',num2str(Temperature),'^3)   + rho^5*(x22/',num2str(Temperature),'^2+x23/',num2str(Temperature),'^4)   + rho^7*(x24/',num2str(Temperature),'^2+x25/',num2str(Temperature),'^3)   + rho^9*(x26/',num2str(Temperature),'^2+x27/',num2str(Temperature),'^4)   + rho^11*(x28/',num2str(Temperature),'^2+x29/',num2str(Temperature),'^3)   + rho^13*(x30/',num2str(Temperature),'^2+x31/',num2str(Temperature),'^3+x32/',num2str(Temperature),'^4) ) * exp(-3.0*rho^2)'), 'independent', {'rho'}, 'dependent', 'P' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.MaxFunEvals = 1000;
%     opts.MaxIter = 10000;
%     opts.Robust = 'Bisquare';

    % Fit model to data.
    [fitresult{indexT}, gof{indexT}] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    hbar = errorbar( xData, yData, Perror(indexT,:), '.', ...
            'LineStyle','none', 'LineWidth', 1.0, ...
            'DisplayName', strcat('$T=',num2str(Temperature,'%.2f'),'$'), ...
            'Parent', axes1 );
    h = plot( fitresult{indexT}, [0; xData], [0; yData] );
%     h = fplot( f(rho), [0 max(DensityRange)*10/9] );
    set( h, 'Parent', axes1, 'Color', get(hbar,'Color'), 'LineWidth', 1.0 );
    
    % Include the errorbar in the legend
    legendVector(indexT) = hbar;
    
%    % Find the roots of the fit
%     fitCoefficients = coeffvalues(fitresult{indexT});
%     f(x) = fitCoefficients(1)*x^5 + fitCoefficients(2)*x^4 + ...
%            fitCoefficients(3)*x^3 + fitCoefficients(4)*x^2 + ...
%            fitCoefficients(5)*x - Psat(indexT);
%     df = diff(f,x);
%     f(rho) = rho.*Temperature+rho.^2.*(fitresult{indexT}.x1.*Temperature+fitresult{indexT}.x2.*sqrt(Temperature)+fitresult{indexT}.x3+fitresult{indexT}.x4./Temperature+fitresult{indexT}.x5./Temperature.^2) ...
%              + rho.^3.*(fitresult{indexT}.x6.*Temperature+fitresult{indexT}.x7+fitresult{indexT}.x8./Temperature+fitresult{indexT}.x9./Temperature.^2) + rho.^4.*(fitresult{indexT}.x10.*Temperature+fitresult{indexT}.x11+fitresult{indexT}.x12./Temperature) ...
%              + rho.^5.*(fitresult{indexT}.x13) + rho.^6.*(fitresult{indexT}.x14./Temperature+fitresult{indexT}.x15./Temperature.^2) + rho.^7.*(fitresult{indexT}.x16./Temperature) ...
%              + rho.^8.*(fitresult{indexT}.x17./Temperature+fitresult{indexT}.x18./Temperature.^2) + rho.^9.*(fitresult{indexT}.x19./Temperature.^2) ...
%              + ( rho.^3.*(fitresult{indexT}.x20./Temperature.^2+fitresult{indexT}.x21./Temperature.^3) ...
%                + rho.^5.*(fitresult{indexT}.x22./Temperature.^2+fitresult{indexT}.x23./Temperature.^4)   + rho.^7.*(fitresult{indexT}.x24./Temperature.^2+fitresult{indexT}.x25./Temperature.^3) ...
%                + rho.^9.*(fitresult{indexT}.x26./Temperature.^2+fitresult{indexT}.x27./Temperature.^4)   + rho.^11.*(fitresult{indexT}.x28./Temperature.^2+fitresult{indexT}.x29./Temperature.^3) ...
%                + rho.^13.*(fitresult{indexT}.x30./Temperature.^2+fitresult{indexT}.x31./Temperature.^3+fitresult{indexT}.x32./Temperature.^4) ) .* exp(-3.0.*rho.^2) - Psat(indexT);
%     df = diff(f,rho);
%    fitZeros{indexT} = vpasolve(f,[min(DensityRange)*1/10 max(DensityRange)*10/9]);

    fitZeros{indexT} = vpasolve(f-Psat(indexT),rho,[0 max(DensityRange)],'random',true);
    if ~isempty(fitZeros{indexT})
        fitZeros{indexT} = sort( [ fitZeros{indexT}(1) vpasolve(f-Psat(indexT),rho,[0 max(DensityRange)],'random',true) ] );
        if length(fitZeros{indexT}) == 2
            if fitZeros{indexT}(1) == fitZeros{indexT}(2) || ( fitZeros{indexT}(2) - fitZeros{indexT}(1) ) < 0.1
                fitZeros{indexT}(2) = [];
            end
        end
    end
    
%     sols=[];
%     for i=1:2
%         sols = [sols double(vpasolve(f-Psat(indexT),rho,[0 max(DensityRange)/2],'random',true))];
%         sols = [sols double(vpasolve(f-Psat(indexT),rho,[max(DensityRange)/2 max(DensityRange)*10/9],'random',true))];
%     end
%     fitZeros{indexT} = unique(sols);

end

    % Create xlabel
    xlabel('$\rho$','FontSize',22,'Interpreter','latex');

    % Create title
    title('$\langle P \rangle$ vs $\rho$ for different isotherms','FontSize',24,...
        'Interpreter','latex');

    % Create ylabel
    ylabel({'$\langle P \rangle$'},'FontSize',22,'Interpreter','latex');

    % Activate box & grid
    box(axes1,'on');
    grid(axes1,'on');

    % Create legend
    legend1 = legend( legendVector );
    set(legend1, 'Location', 'BestOutside', 'Interpreter','latex','FontSize',18);
    

%% Plot the Maxwell construction results
% Get the points
TvsRho = [];
TvsRhoErrors = [];
for indexT = 1:length(TemperatureRange)
    temp = TemperatureRange(indexT);
    if length(fitZeros{indexT}) == 2
        zero1 = double(fitZeros{indexT}(1));
        zero2 = double(fitZeros{indexT}(2));
    elseif length(fitZeros{indexT}) == 1
         zero1 = double(fitZeros{indexT}(1));
         zero2 = 0;
    end
    if ( length(fitZeros{indexT}) == 2 || length(fitZeros{indexT}) == 1 ) && temp < 1.35
        zeroerror1 = 0.0;
        zeroerror2 = 0.0;
%         % Use the poly5 fit
%         fitCoefficients = coeffvalues(fitresult{indexT});
%         f(x) = fitCoefficients(1)*x^5 + fitCoefficients(2)*x^4 + ...
%                fitCoefficients(3)*x^3 + fitCoefficients(4)*x^2 + ...
%                fitCoefficients(5)*x - Psat(indexT);
%         df = diff(f,x);
% 
%         zeroPerror = diff(predint(fitresult{indexT},[min([zero1 zero2]) max([zero1 zero2])],0.95,'functional','off')')./2;
%         zerodP1 = double(df(zero1));
%         zerodP2 = double(df(zero2));
%         zeroerror1 = zeroPerror(1)/zerodP1;
%         zeroerror2 = zeroPerror(2)/zerodP2;
%         zeroRhoSI1 = linspace(min(DensityRange)*1/10,...
%                   min([zero1,zero2])+abs(zero1+zero2)/2,100);
%         zeroRhoSI2 = linspace(max([zero1,zero2])-abs(zero1+zero2)/2,...
%                   max(DensityRange)*12/9,100);
%         zeroPCI1 = predint(fitresult{indexT}, zeroRhoSI1,...
%                    0.95,'functional','off')./2;
%         zeroPCI2 = predint(fitresult{indexT}, zeroRhoSI2,...
%                   0.95,'functional','off')./2;
%         zeroerror1 = zeroerror1 + findRhoError(zero1, zeroRhoSI1, zeroPCI1);
%         zeroerror2 = zeroerror2 + findRhoError(zero2, zeroRhoSI2, zeroPCI2);
        
        zeroerror1 = 0.015;
        zeroerror2 = 0.015;

        TvsRho = [ TvsRho; zero1 temp];
        TvsRhoErrors = [TvsRhoErrors; zeroerror1];
        if length(fitZeros{indexT}) == 2
             TvsRho = [TvsRho; zero2 temp ];
             TvsRhoErrors = [TvsRhoErrors; zeroerror2];
        end
    end
end

% Create figure
figure2 = figure('PaperOrientation','landscape','PaperType','A3');

% Create axes
axes2 = axes('Parent',figure2);
hold(axes2,'on');

% Plot Johnson's results
johnsonLine = plot([JohnsonData(:,2); JohnsonData(:,3)],[JohnsonData(:,1); JohnsonData(:,1)],'DisplayName','Johnson et al.',...
    'Parent',axes2,...
    'LineStyle','none',...
    'LineWidth',1.5,...
    'MarkerSize',4,...
    'MarkerFaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'MarkerEdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'Marker','o');

% Create plot
plot2 = errorbar(TvsRho(:,1),TvsRho(:,2),TvsRhoErrors,'horizontal',...
    'DisplayName','Me','MarkerSize',4,...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'Marker','square',...
    'LineStyle','none',...
    'LineWidth',1.5,...
    'Color',[0 0.447058826684952 0.74117648601532]);

% Get xdata from plot
xdata2 = get(plot2, 'xdata');
% Get ydata from plot
ydata2 = get(plot2, 'ydata');
% Make sure data are column vectors
xdata2 = xdata2(:);
ydata2 = ydata2(:);

% Remove NaN values and warn
nanMask1 = isnan(xdata2(:)) | isnan(ydata2(:));
if any(nanMask1)
    warning('GeneratedCode:IgnoringNaNs', ...
        'Data points with NaN coordinates will be ignored.');
    xdata2(nanMask1) = [];
    ydata2(nanMask1) = [];
end

% Find x values for plotting the fit based on xlim
axesLimits2 = xlim(axes2);
xplot1 = linspace(axesLimits2(1), axesLimits2(2), 1000);

% Set up fittype and options.
fittypePhase = fittype( 'poly3' );
optsPhase = fitoptions( 'Method', 'LinearLeastSquares' );
optsPhase.Robust = 'Bisquare';
%optsPhase.Weights = 1./TvsRhoErrors.^2;

% Do the fit
fitResults1 = fit( xdata2, ydata2, fittypePhase, optsPhase );

% Find coefficients for polynomial (order = 4)
%fitResults1 = fit(xdata2,ydata2,'poly4');

% Evaluate polynomial
yplot1 = feval(fitResults1,xplot1);

% Find the maximum
[MaxValue, MaxIdx] = max(yplot1);
MaxRho = xplot1(MaxIdx(1));

% Find the errorbars of the maximum points
MaxTPredInt = predint(fitResults1,MaxValue(1),0.95,'functional','off');
MaxValueError = abs(diff(MaxTPredInt'))/2;
RhoSearchInt = linspace(min(DensityRange), max(DensityRange), 1000);
CISearchInt = predint(fitResults1,RhoSearchInt,0.95,'functional','off');
MaxRhoError = findEquals(CISearchInt(:,2), RhoSearchInt, feval(fitResults1,RhoSearchInt));

% Put everything in vectors
MaxValue = [MaxValue(1) MaxRho];
MaxValueError = [MaxValueError MaxRhoError];

% % Plot the fit
% fitLine1 = plot(xplot1,yplot1,'DisplayName','3rd degree poly fit',...
%     'Tag','3rd degree',...
%     'Parent',axes2,...
%     'LineWidth',1.5,...
%     'Color',[0 0.447058826684952 0.74117648601532]);
% 
% % Set new line in proper position
% setLineOrder(axes2,fitLine1,plot2);

% Create xlabel
xlabel('$\rho$','FontSize',22,'Interpreter','latex');

% Create title
title('Coexistence curve of the LJ fluid','FontSize',24,...
    'Interpreter','latex');

% Create ylabel
ylabel('$T$','FontSize',22,'Interpreter','latex');

box(axes2,'on');
% Set the remaining axes properties
set(axes2,'FontSize',18,'XGrid','on','YGrid','on');

% % Do a crazy fit and plot it
%crazyfitresult = createCrazyFit(TemperatureRange, DensityRange, Pmean-Psat');
% Plot Johnson
% clear clazyfitresult;
% crazyfitresult = JohnsonParams;
%  fimplicit( @(rho,T) rho.*T+rho.^2.*(crazyfitresult.x1.*T+crazyfitresult.x2.*sqrt(T)+crazyfitresult.x3+crazyfitresult.x4./T+crazyfitresult.x5./T.^2) ...
%               + rho.^3.*(crazyfitresult.x6.*T+crazyfitresult.x7+crazyfitresult.x8./T+crazyfitresult.x9./T.^2) + rho.^4.*(crazyfitresult.x10.*T+crazyfitresult.x11+crazyfitresult.x12./T) ...
%               + rho.^5.*(crazyfitresult.x13) + rho.^6.*(crazyfitresult.x14./T+crazyfitresult.x15./T.^2) + rho.^7.*(crazyfitresult.x16./T) ...
%               + rho.^8.*(crazyfitresult.x17./T+crazyfitresult.x18./T.^2) + rho.^9.*(crazyfitresult.x19./T.^2) ...
%               + ( rho.^3.*(crazyfitresult.x20./T.^2+crazyfitresult.x21./T.^3) ...
%                 + rho.^5.*(crazyfitresult.x22./T.^2+crazyfitresult.x23./T.^4)   + rho.^7.*(crazyfitresult.x24./T.^2+crazyfitresult.x25./T.^3) ...
%                 + rho.^9.*(crazyfitresult.x26./T.^2+crazyfitresult.x27./T.^4)   + rho.^11.*(crazyfitresult.x28./T.^2+crazyfitresult.x29./T.^3) ...
%                 + rho.^13.*(crazyfitresult.x30./T.^2+crazyfitresult.x31./T.^3+crazyfitresult.x32./T.^4) ) .* exp(-crazyfitresult.g.*rho.^2) ,...
%             [min(DensityRange) max(DensityRange)*10/9 min(TemperatureRange)*55/60 max(TemperatureRange)],...
%             'DisplayName','MBWR fit',...
%             'Parent',axes2,...
%             'LineWidth',1.5);
       
% Create legend
legend(axes2,'show');


%-------------------------------------------------------------------------%
function setLineOrder(axesh1, newLine1, associatedLine1)
%SETLINEORDER(AXESH1,NEWLINE1,ASSOCIATEDLINE1)
%  Set line order
%  AXESH1:  axes
%  NEWLINE1:  new line
%  ASSOCIATEDLINE1:  associated line

% Get the axes children
hChildren = get(axesh1,'Children');
% Remove the new line
hChildren(hChildren==newLine1) = [];
% Get the index to the associatedLine
lineIndex = find(hChildren==associatedLine1);
% Reorder lines so the new line appears with associated data
hNewChildren = [hChildren(1:lineIndex-1);newLine1;hChildren(lineIndex:end)];
% Set the children:
set(axesh1,'Children',hNewChildren);


%-------------------------------------------------------------------------%
function RhoError = findRhoError(rhoZero, rhoInterval, PCInterval)

% Determine approximately the position of the zero inside the given
% interval
zeroIdx = round((rhoZero-min(rhoInterval))/(max(rhoInterval)-min(rhoInterval))*length(rhoInterval));

% Find the right bound
diffPrevUp = abs(PCInterval(zeroIdx,2));
diffPrevLo = abs(PCInterval(zeroIdx,1));
for index = zeroIdx+1:length(rhoInterval)
    diffUp = abs(PCInterval(index,2));
    diffLo = abs(PCInterval(index,1));
    if diffUp > diffPrevUp || diffLo > diffPrevLo
        rightBound = rhoInterval(index-1);
        break
    end
    diffPrevUp = diffUp;
    diffPrevLo = diffLo;
end

% Find the left bound
diffPrevUp = abs(PCInterval(zeroIdx,2));
diffPrevLo = abs(PCInterval(zeroIdx,1));
for index = zeroIdx-1:-1:1
    diffUp = abs(PCInterval(index,2));
    diffLo = abs(PCInterval(index,1));
    if diffUp > diffPrevUp || diffLo > diffPrevLo
        leftBound = rhoInterval(index+1);
        break
    end
    diffPrevUp = diffUp;
    diffPrevLo = diffLo;
end

RhoError = (rightBound-leftBound)/2;


%-------------------------------------------------------------------------%
function RhoError = findEquals(upperCI, RhoInt, TInt)

[MaxT, MaxIdx] = max(TInt);

diffPrev = abs(upperCI(1)-MaxT);
for index = 2:MaxIdx
    diff = abs(upperCI(index)-MaxT);
    if diff > diffPrev
        leftBound = RhoInt(index-1);
        break
    end
    diffPrev = diff;
end

diffPrev = abs(upperCI(MaxIdx)-MaxT);
for index = MaxIdx+1:length(RhoInt)
    diff = abs(upperCI(index)-MaxT);
    if diff > diffPrev
        rightBound = RhoInt(index-1);
        break
    end
    diffPrev = diff;
end

RhoError = (rightBound-leftBound)/2;
    
    

