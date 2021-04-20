function [cr, pvalue] = CC_inference(varargin)
%
% CC analysis (select for surrogates or th-ing)
%
%
% Inputs:
%varargin{1} is the Original time - series. Dimensions: 
%varargin{2} is the network of threshold (obtained from Surrogates) for each link.
%varargin{3} is the dimension of the data: 
%           set it to 2 for concatenated all
%           set it to 3 for subject
%           set it to 4 for recording-subject.
% Outputs: 
%       cr => Partial Correlation network, thresholded.
%       pvalue => p-value that comes from the output of the CC built-in function 
%                 You can ignore this since it is only useful for thresholding and this program uses Surrogates to threshold.
%
A = varargin{1};
thr = varargin{2}; % 2 dimensional matrix of thresholds 

if varargin{3} == 2
    % concatenated
    
    % 1. get the correlations
    [rho, pvalue] = CalculofCC(A,2);
    % 2. threshold
    cr = rho.*(rho >= thr);
    
    
elseif varargin{3} == 3
    % subject (S) network
    A = varargin{1};
    [rho, pvalue] = CalculofCC(A,3);
    % 2. threshold
    cr = zeros(size(rho,1), size(rho,2), size(rho,3));
    for i = 1:size(rho,3)
        tmp = rho(:,:,i)
        cr = tmp.*(tmp >= thr(:,:,i));
    end
       
elseif varargin{3} == 4
    % recording-subject (RS) network
     A = varargin{1};
    [rho, pvalue] = CalculofCC(A,4);
    % 2. threshold
    cr = zeros(size(rho,1), size(rho,2), size(rho,3), size(rho,4));
    for i = 1:size(rho,3)
        for j = 1:size(rho,4)
            tmp = rho(:,:,i,j);
            cr = tmp.*(tmp >= thr(:,:,i,j));
        end
    end
       
end


end
