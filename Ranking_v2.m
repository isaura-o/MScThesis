function [Y3, bL, xcor, ycor, storC] = Ranking_v2(varargin)

A = abs(varargin{1}); % concatenate 1 (2dim mat)
B = abs(varargin{2}); % concatenate 2 (2dim mat)
A1 = abs(varargin{3}); % Intra-links network for method 1 (2dim mat)
B1 = abs(varargin{4}); % Intra-links network for method 2 (2dim mat)


% concatenate
tmpC = A(~tril(ones(size(A))));
% concatenate
tmp = B(~tril(ones(size(B))));

tmp1 = A1(~tril(ones(size(A1))));
tmp2 = B1(~tril(ones(size(B1))));

% sort the values (put another 0 afterwards...)
[valC,numC] = sort(tmpC, 'descend');
storC = numC;

[val,num] = sort(tmp, 'descend');
stor = num;
% the values have labels 1 to 36 (66 or lenght of the vector tmp) ( rank number)
lab = (1:length(tmpC))';

% matrix of the labels and the position of vector 2 = [ label, linknum ];
% impose that al the valPCXX == 0; must have the label or linknum == 0
lab( valC == 0 ) = 0;
num(val == 0 ) = 0;

% find the values that match. 
[~, pos] = ismember( numC, num);

X = [lab, pos];
T1 = X(:,1) == 0 & X(:,2) == 0;
Y = X;
Y(T1,:)= [];


% valors off diagonal (non block diag)
[val2,num2] = sort(tmp1, 'descend');
[val2a,num2a] = sort(tmp2, 'descend');

numsnondig = num2(val2>0);
numsnondig2 = num2a(val2a>0);

% find the initial positions on the sorted full matrices
% this gives me the initial positions of the nondiagonal values
a = zeros(1,length(numsnondig));
for i = 1:length(numsnondig)
    a(i) = find(numC == numsnondig(i));
end
b = zeros(1,length(numsnondig2));
for i = 1:length(numsnondig2)
    b(i) = find(num == numsnondig2(i));
end

% find if the non diagonal values are a zero or non zero in the final
% positions

% for the MECMI that is the "base" so ranks 1 to x taking into account that
% 1 is sorted the strongest.
% for this we search which is the position of the value in the 1 to x nums
% and set them 1 for non diagonal value
x = sum((lab == a),2);

% find the values of the vector pos that get the initial values
ndx = pos == 0;
vec(ndx) = 0; vec(~ndx) = val(pos(~ndx));
vec2 = vec';

% now for the other values we need to take into account that they are
% sorted respect the base matrix, so we need to find the positions in the
% ranking
[~,x2] = ismember(numC, num(b));
x2a = x2>0;
X2 = [lab, pos, x2a, x, valC, vec2];
T12 = X2(:,1) == 0 & X2(:,2) == 0 & X2(:,3) == 0 & X2(:,4) == 0; 
Y2 = X2;
Y2(T12,:)= [];
Y3 = [Y2(:,1) Y2(:,2) (Y2(:,3)|Y2(:,4)) Y2(:,5) Y2(:,6)];


% now give the values in Y3 where the column 1 is a number but column 2 is
% a zero
blackList = [Y3(Y3(:,2) == 0 ), Y3( Y3(:,2) == 0 , 4)];
% tell me the coordinates in the network of this ones
xcor = zeros(2,size(blackList,1));
ycor = zeros(2,size(blackList,1));
 for i = 1:size(blackList,1)
     [ycor(:,i),xcor(:,i)] = find (blackList(i,2) == A);
 end
 
 bL = blackList;

end
