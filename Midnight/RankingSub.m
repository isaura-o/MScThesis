function [rancir, indxStrZ, maxRankOnY, maxRankOnY5] = RankingSub(varargin)
% rank the subjects, taking x strongest links.

X = varargin{1};
Y = varargin{2};
X1 = varargin{3};
Y1 = varargin{4};
z = varargin{5};

%1. use ranking function.
[Y3, ~, ~, ~, sortC] = Ranking_v2(X,Y,X1,Y1); %Y3 has the indexes needed for indexing

%2. find strongest links: take only the zth first values
rancir = Y3(1:z,:);
indxStrZ = sortC(1:z,:);
maxRankOnY = max(Y3(:,2));
maxRankOnY5 = max(Y3(1:z,2));
