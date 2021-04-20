function [Y, stor] = Ranking(varargin)

A = abs(varargin{1}); % concatenate (2dim mat)
B = abs(varargin{2}); % most probable network (2dim mat)

% concatenate
tmpC = A(~tril(ones(size(A))));
% most probable network
tmp = B(~tril(ones(size(B))));

% sort the values (put another 0 afterwards...)
[valC,numC] = sort(tmpC, 'descend');
stor = numC;

[val,num] = sort(tmp, 'descend');
% the values have labels 1 to 36 (66 or lenght of the vector tmp) ( rank number)
lab = (1:length(tmpC))';

% matrix of the labels and the position of vector 2 = [ label, linknum ];
% impose that al the valPCXX == 0; must have the label or linknum == 0
lab( valC == 0 ) = 0;
num(val == 0 ) = 0;

% find the values that match. 
[~, pos] = ismember( numC, num);


% find the values of the vector pos that get the initial values
ndx = pos == 0;
vec(ndx) = 0; vec(~ndx) = val(pos(~ndx));
vec2 = vec';

X = [lab, pos, valC, vec2];
T1 = X(:,1) == 0 & X(:,2) == 0 & X(:,3) == 0 & X(:,4) == 0;
Y = X;
Y(T1,:)= [];


end
