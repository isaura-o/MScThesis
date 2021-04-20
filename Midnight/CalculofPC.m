function [prho, pvalue] = CalculofPC(varargin)
%
% Calculation of cross correlation. Works for data and surrogate data.
%
% input:
% X = time series
% 
% dim = dimension of the time series.
%       2 =>ENTER: nodes x timepoints.
%       3 =>ENTER: nodes x timepoints x subjects.
%       4 =>ENTER  nodes x timepoints x recordings x subjects.
% 
%
% NOTE: uncomment "progressbar" function if you want to control the time it lasts.
%
X = varargin{1};
dim = varargin{2};

if (dim == 2)
    Y = reshape(X, size(X,1),[]);
    Y = Y - mean(Y,2);
    [prho,pvalue] = partialcorr(Y');  prho = abs(prho - eye(size(Y,1))); 

elseif (dim == 3)
    X = permute(X, [1 4 2 3]); 
    X1c = reshape(X, size(X,1), size(X,2), []);
    X1c = X1c - mean(X1c, 3);
    Y = permute(X1c, [1 3 2]);
    
    prho = zeros(size(Y,1), size(Y,1), size(Y,3));
    pvalue = zeros(size(Y,1), size(Y,1), size(Y,3));
%    progressbar
    for i = 1:size(Y,3)
       X = Y(:,:,i);
       [tmp, pvalue(:,:,i)] = partialcorr(X');
       prho(:,:,i) = abs(tmp - eye(size(tmp,1)));     
%       progressbar(i/size(Y,3))
    end
        
elseif (dim == 4)
    Y = X;
    prho = zeros(size(Y,1), size(Y,1), size(Y,3), size(Y,4));
    pvalue = zeros(size(Y,1), size(Y,1), size(Y,3), size(Y,4));
%   progressbar
    for i = 1:size(Y,4)
        for j = 1:size(Y,3)
            X = Y(:,:,j,i);
            [tmp, pvalue(:,:,j,i)] = partialcorr(X');
            prho(:,:,j,i) = abs(tmp - eye(size(Y,1)));
        end
%       progressbar(i/size(Y,4))
    end
      
elseif (dim == 5) 
    % Only useful for Surrogates and 5 min slices 5-dim arrays, don't use otherwise.
    Y = X;
    prho = zeros(size(Y,1), size(Y,1), size(Y,3), size(Y,4), size(Y,5));
%   progressbar
    for i = 1:size(Y,3)
        for j = 1:size(Y,4)
            for k = 1:size(Y,5)
                X = Y(:,:,i,j,k);
                [tmp, ~] = partialcorr(X');
                prho(:,:,i,j,k) = (abs(tmp) - eye(size(Y,1)));
            end
        end
%       progressbar(i/size(Y,3))
    end
end


end


