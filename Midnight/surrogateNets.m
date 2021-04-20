function [rhoTh, prhoTh] = surrogateNets(varargin)
%
%
% Surrogate network, uses the function surrogates.
%
% varargin{1} = Time-series. (dimensions: ROIxTSxRecordingxSubject)
% varargin{2} = dimension of the surrogates:
%               varargin{2} == 2 => ROIxTS (concatenated)
%               varargin{2} == 3 => ROIxTSxSubject (Recordings concatenated)
%               varargin{2} == 4 => ROIxTSxRecordingsxSubject
%
% confidence interval, you can set it at 0.95 in main.m)
Th = varargin{2};

% Load the Time-series
X = varargin{1};

if varargin{3} == 2
    %concatenate
    X1 = reshape(X, size(X,1),[]);
    
    rho = zeros(size(X1,1), size(X1,1), 1000);
    prho = zeros(size(X1,1), size(X1,1), 1000);
    for k = 1:1000
            % surrogate time-series
            Yc = surrogates(X1);
            % surrogate CC networks
            rho(:,:,k) = CalculofCC(Yc,2);
            % surrogate PC networks
            prho(:,:,k) = CalculofPC(Yc,2);
    end
    clear Yc X1 XS
    
    %get the thresholds for each link in each SRnetwork
    rhoTh = thres9095_v4(rho,Th,1);
    prhoTh = thres9095_v4(prho,Th,1);

elseif varargin{3} == 3
     %Subject networks
     %This step is not really necessary but preferable since it may reshape badly if we directly 
     %reshape dims 2 and 3 together.
     X = permute(X, [1 4 2 3]); 
     
     X1c = reshape(X, size(X,1), size(X,2), []);
     XS = permute(X1c, [1 3 2]);

     rho = zeros(size(XS,1), size(XS,1), 1000, size(XS,3));
     prho = zeros(size(XS,1), size(XS,1), 1000, size(XS,3));
     for j = 1:size(XS,3)
          for k = 1:1000
              % surrogate time-series
              Y = surrogates(XS(:,:,j));
              % surrogate CC networks
              rho(:,:,k,j) = CalculofCC(Y,2);
              % surrogate PC networks
              prho(:,:,k,j) = CalculofPC(Y,2);
          end
     end
    clear Yc X1c
    %get the thresholds for each link in each SRnetwork
    rhoTh = thres9095_v4(rho,Th,2);
    prhoTh = thres9095_v4(prho,Th,2);
    
%    clear rho prho
    
elseif varargin{3} == 4
%   Surrogates for Subject and Recording.

    rho = zeros(size(X,1), size(X,1), 1000, size(X,3), size(X,4));
    prho = zeros(size(X,1), size(X,1), 1000, size(X,3), size(X,4));
    for i = 1:size(X,3)
        for j = 1:size(X,4)
            for k = 1:1000
                Y = surrogates(X(:,:,j,i));
                rho(:,:,k,j,i) = CalculofCC(Y,2);
                prho(:,:,k,j,i) = CalculofPC(Y,2);
            end
        end
    end
    clear Y XS
    
    %get the thresholds for each link in each SRnetwork
    rhoTh = thres9095_v4(rho,Th,3);
    prhoTh = thres9095_v4(prho,Th,3);
    
%    clear rho prho
  
    
end

end