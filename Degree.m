function [degC, degS, aT, bT] = Degree(varargin)
%
% Degree for concatenate and subjects
%
% Distribution of links for subjects.
% 
% set dimensionality in varargin = 2 for obtaining distributions for degree in subjects and
% in recordings ---- if not needed set dim = 1
%
% Inputs:
%        varargin{1} = Concatenate networks
%        varargin{2} = Subject networks
%
% Output:
%       degC = degree for concatenate network
%       degS = degree for Subject networks
%       degRS = degree for Recording-subject networks (you may not want it)
%

% concatenated
degC = sum(sum(triu(logical(varargin{1}),1)));

% subjects
Y = varargin{2};
degS = zeros(1,size(Y,3)); 
for k = 1:size(Y,3)
	degS(k) = sum(sum(triu(logical(Y(:,:,k)),1)));
end

% distribution for all subjects
maxD = max(degS(:)); minD = min(degS(:));
cntrs = minD:maxD;
[aT, bT] = histogram(degS(:), cntrs);

% recording-subject degree 
% Y = varargin{3};
% degRS = zeros(size(Y,3),size(Y,4), size(Y,5)); 
% for i = 1:size(varargin{3},3)
%     for k = 1:size(varargin{3},4)
%         for j = 1:size(varargin{3},5)
%         	degRS(i,k,j) = sum(sum(triu(logical(Y(:,:,i,k,j)),1)));
%         end
%     end
% end    
         
         
    

end