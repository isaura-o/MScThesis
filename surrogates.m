function [Y] = surrogates(X)
% Calls Fourier algorithm to do the surrogate serie: iaaft_loop_1d.m  
%
% Inputs
% Y = surrogate data
%
% Outputs
% X = Time series. 
% 
N = size(X,1);
    
Y = zeros(N,size(X,2));
for i = 1:N
	%1. mean value
    meanValue = mean(X(i,:));
    %1. ordered list of the TS(time series):
    srtX = sort(X(i,:) - meanValue);
    %2. amplitudes of TF of the original TS(time series):
    s2 = abs(ifft(X(i,:) - meanValue)); 
    %Using iaaft_loop_1d:
    fourierCoeff = s2;
    sortedValues = srtX;
    [X1, errAmplitud, errSpc] = iaaft_loop_1d(fourierCoeff, sortedValues);
    
    Y(i,:) = X1;
end
 
end
