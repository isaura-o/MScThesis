function [rho, pvalue] = CalculofCC(varargin)
%
% Calculation of cross correlation. Works for data and surrogate data.
%
% input:
% Y = time series
% 
% dim = dimension of the time series. 
%       2 => nodes x timepoints.
%       3 => nodes x timepoints x subjects (or num of surrogates).
%       4 => nodes x timepoints x recordings x subjects(or num of surrogates).
%       5 => nodes x timepoints x recordings x subjects x num of surrogates).
% 
% Outupts:
% rho = CC network.
% pvalue
%
X = varargin{1};
dim = varargin{2};

if (dim == 2)
    Y = reshape(X, size(X,1),[]);
    Y = Y - mean(Y,2);
    [rho,pvalue] = corrcoef(Y'); rho = abs(rho - eye(size(Y,1)));
        
elseif (dim == 3)
    X = permute(X, [1 4 2 3]); 
    X1c = reshape(X, size(X,1), size(X,2), []);
    X1c = X1c - mean(X1c, 3);
    Y = permute(X1c, [1 3 2]);
    
    rho = zeros(size(Y,1), size(Y,1), size(Y,3));
    pvalue = zeros(size(Y,1), size(Y,1), size(Y,3));
    for i = 1:size(Y,3)
        [tmp, pvalue(:,:,i)] = corrcoef(Y(:,:,i)'); rho(:,:,i) = abs(tmp - eye(size(tmp,1)));
    end
        
elseif (dim == 4)
    Y = X;
    rho = zeros(size(Y,1), size(Y,1), size(Y,3), size(Y,4));
    pvalue = zeros(size(Y,1), size(Y,1), size(Y,3), size(Y,4));
    for i = 1:size(Y,4)
        for j = 1:size(Y,3)
            [tmp, pvalue(:,:,j,i)] = corrcoef(Y(:,:,j,i)');
            rho(:,:,j,i) = abs(tmp - eye(size(Y,1)));
        end
    end
      
elseif (dim == 5)
    %Only useful for surrogates don't give pvalue
    Y = X;
    rho = zeros(size(Y,1), size(Y,1), size(Y,3), size(Y,4), size(Y,5));
    for i = 1:size(Y,3)
        for j = 1:size(Y,4)
            for k = 1:size(Y,5)
                [tmp,~] = corrcoef(Y(:,:,i,j,k)');
                rho(:,:,j,i) = abs(tmp - eye(size(Y,1)));
            end
        end
    end
    
end

end

