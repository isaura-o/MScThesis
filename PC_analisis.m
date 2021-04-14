function [pr, pvalue] = PC_analisis(varargin)
%
% PC analysis (select for surrogates or th-ing) (Same as CC_analisisBD_v2)
%
% 1. Concatenate (A = concatenated data set) 
% 2. Subjects (B = subject data set), (rhos = non concatenated surrogate networks)
% 3. Threshold (needed input of matrix of thresholds: thr)
%
% Inputs:
%varargin{1} is the Original time - series.
%varargin{2} is the network of threshold (obtained from Surrogates) for each link.
%varargin{3} is the dimension of the data: 
%           set it to 2 for concatenated all
%           set it to 3 for other
% Outputs: 
%       pr => Partial Correlation network, thresholded.
%       pvalue => p-value that comes from the output of the PC built-in function 
%                 You can ignore this since it is only useful for thresholding and this program uses Surrogates to threshold.
%
A = varargin{1};
thr = varargin{2}; % 2 dimensional matrix of thresholds 

if varargin{3} == 2
    % concatenated
    
    % 1. get the correlations
    [prho, pvalue] = CalculofPC(A,2);
    % 2. threshold
    pr = prho.*(prho >= thr);
    
    
elseif varargin{3} == 3
    % subject (S) network
    A = varargin{1};
    [prho, pvalue] = CalculofPC(A,3);
    % 2. threshold
    pr = zeros(size(prho,1), size(prho,2), size(prho,3));
    for i = 1:size(prho,3)
        tmp = prho(:,:,i)
        pr = tmp.*(tmp >= thr(:,:,i));
    end
       
elseif varargin{3} == 4
    % recording-subject (RS) network
     A = varargin{1};
    [prho, pvalue] = CalculofPC(A,4);
    % 2. threshold
    pr = zeros(size(prho,1), size(prho,2), size(prho,3), size(prho,4));
    for i = 1:size(prho,3)
        for j = 1:size(prho,4)
            tmp = prho(:,:,i,j);
            pr = tmp.*(tmp >= thr(:,:,i,j));
        end
    end
       
end


end
