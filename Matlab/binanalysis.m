function Verror = binanalysis( data, doplot )
%BINANALYSIS gives the error on the potential via the binning technique.
%   Detailed explanation goes here

    % Define the graining of the Lbin scale
%    Nbin_powers_interval = 10.^(2:log10(length(data)));
%    Nbin_powers_interval = 10.^(floor(log10(length(data)))-3:floor(log10(length(data)))-1);
    Nbin_powers_interval = 10.^(1.5:floor(log10(length(data)))-3);
    Nbin_internal_interval = 1:0.01:9;
    Nbin_interval = floor(kron(Nbin_powers_interval, Nbin_internal_interval));
    Nbin_interval = [30:1:1000 1050:50:10000 10500:500:20000] ;
    Length_interval = floor(length(data)./Nbin_interval);

    % Do a preallocation of the vector for speed
%    iteration = 1;
    x_errors = zeros(1,length(Nbin_interval));

%     % Calculate the Vmean and its std for every step
%     for Nbin = Nbin_interval
%         data_binned = data(1:Nbin*Length_interval(iteration));
%         data_binned = reshape(data_binned,[Length_interval(iteration),Nbin]);
%         x_binned_mean = mean(data_binned);
%         bin_std = std(x_binned_mean);
%         x_errors(iteration) = bin_std/sqrt(Nbin);
%         iteration = iteration + 1;
%     end
    
    for NumBinIdx = 1:length(Nbin_interval)
        fprintf( 'Binning technique: doing the calculations for %d bins.\n', Nbin_interval(NumBinIdx));
        x_binned_mean = zeros(1,Nbin_interval(NumBinIdx));
        for BinIdx = 1:Nbin_interval(NumBinIdx)
            x_binned_mean(BinIdx) = sum(data(1+Length_interval(NumBinIdx)*(BinIdx-1):Length_interval(NumBinIdx)*BinIdx))/Length_interval(NumBinIdx);
        end
        bin_std = std(x_binned_mean);
        x_errors(NumBinIdx) = bin_std/sqrt(Nbin_interval(NumBinIdx));
    end

    % Calculate the error as the max of the errors
    Verror = max(x_errors);

    if ( strcmp(doplot, 'Yes') )
        % Create figure
        BinFigure = figure('PaperOrientation','landscape','PaperType','A3');

        % Create axes
        BinAxes = axes('Parent',BinFigure);
        hold(BinAxes,'on');

        % Create plot
        plot(Length_interval, x_errors,'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
            'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
            'MarkerSize',4,...
            'Marker','square',...
            'LineWidth',2,...
            'LineStyle','none');

        % Create xlabel
        xlabel('Bin length','FontSize',24,'Interpreter','latex');

        % Create title
        title('$\sigma_V$ vs bin length','FontSize',24,'Interpreter','latex');

        % Create ylabel
        ylabel('$\sigma_V$ [$\varepsilon$]','FontSize',24,'Interpreter','latex');

        box(BinAxes,'on');
        % Set the remaining axes properties
        set(BinAxes,'FontSize',18);%,'YMinorTick','on','YScale','log');
    end
    

end

