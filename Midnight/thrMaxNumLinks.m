function X = thrMaxNumLinks(varargin)
%
% Give the thresholded network when we do specify the maximum strongest
% links we want.
% example: z = 25 will give us a network with the 25 strongest links.
%
% input:
% A = networks to threshold
% z = max num of links to be on network (threshold)
%
% output:
% X = network threshold
%
%

A = varargin{1};
z = varargin{2};

% make network triangular (so don't get repetitions)
A2 =  triu(A,1);
% sort values
[val, ~] = sort(A2(:), 'descend');
thr = val(z);
% purge values
B = (A2 + A2'); 
X = B.*(B >= thr);

end