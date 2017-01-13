function [ data, gr, Vmean, Verror, Cv, Cverror, Pmean, Perror ] = SingleDataAnalysis( outputfile, groutputfile, temperature , Np)
%SingleDataAnalysis Analyze a single dataset produced by LJMC
%   Imports the data, plots the distribution of the potential values,
%   calculates the error on the mean value of the potential via the
%   binning technique, and calculates the specific heat at constant volume.
%
%
% Example:
%   [data gr Vmean Verror Cv Cverror Pmean Perror] = SingleDataAnalysis('outfile', 'g(r)_outfile', 2.0, 100);


    % Import the data
    data = importdata(outputfile);
    %load workspace.mat
    gr = importdata(groutputfile);
    
    % Save the workspace, because it could crash later
    save workspace
    
    % Plot the distribution of the potential values
    [~, Vmean, ~] = plotdistr(data(:,1)/Np);
    
    % Plot the behavior of the mean value against the MCTime
    plotmctime(data(:,1)/Np);
    
    % Estimate the error via the binning analysis
    Verror = binanalysis(data(:,1)/Np, 'Yes');
    
    % Calculate Cv and its error
    [Cv, Cverror] = specificheat(data(:,1), temperature, Np);
    
    % Calculate P and its error
    Pmean = mean(data(:,2));
    Perror = binanalysis(data(:,2), 'No');
    
    % Plot the autocorrelation; DEACTIVATE for large sets
    %plotautocorrelation(data(:,1)/Np);
    
    % Plot the radial distribution function
    plotrdf(gr);
    
    % Print out some infos
    fprintf( '\n V = %.4f +/- %.4f \n Cv = %.4f +/- %.4f \n\n', Vmean, Verror, Cv, Cverror);

end

