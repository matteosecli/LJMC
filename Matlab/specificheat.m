function [ Cv, deltaCv ] = specificheat( data, temp, np )
%SPECIFICHEAT Summary of this function goes here
%   Detailed explanation goes here

    % Calculate the specific heat
    Cv = (mean(data.*data) - mean(data)^2)/(temp^2)/np;
    
    %% Estimate its error by dividing in bins.
    % The BinWidth should be around the point at which the autocorrelation 
    % goes to zero. It seems to be 10E4 in my simulation.
    BinWidth = 10000;
    Nbins = length(data)/BinWidth;
    
    % Do a preallocation of the vector for speed
    iteration = 1;
    Cv_vector = zeros(1,Nbins);
    
    % Calvulate the Cv for each bin
    data_binned = reshape(data,[BinWidth,Nbins]);
    for bin = 1:1:Nbins
        Cv_vector(iteration) = (mean(data_binned(:,bin).*data_binned(:,bin)) - mean(data_binned(:,bin))^2)/(temp^2);
        iteration = iteration + 1;
    end

    % Calculate the error on Cv as the standard deviation of the Cv's in
    % the bin, renormalized by the number of bins.
    deltaCv = std(Cv_vector)/sqrt(Nbins)/np;
    
end

