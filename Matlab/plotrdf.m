function plotrdf(gr)
%PLOTRDF(GR)
%  GR: matrix that contains r as the first column and g(r) as the second
%  one.


% Create figure
figure1 = figure('PaperOrientation','landscape','PaperType','A3');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(gr(:,1),[gr(:,2), ones(length(gr(:,1)),1)],'Parent',axes1);
set(plot1(1),'MarkerSize',4,'Marker','square','LineWidth',2.5);
set(plot1(2),'LineWidth',2.5,'LineStyle','--',...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625]);

% Create xlabel
xlabel('$r$ [$\sigma$]','FontSize',24,'Interpreter','latex');

% Create title
title('Radial distribution function','FontSize',24,'Interpreter','latex');

% Create ylabel
ylabel('$g(r)$','FontSize',24,'Interpreter','latex');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 gr(end,1)]);

box(axes1,'on');

% Set the remaining axes properties
set(axes1,'FontSize',18);
