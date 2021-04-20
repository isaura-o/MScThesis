function [a1, b1, val] = locNodes(varargin)
% 
% This function finds the coordinates in the matrix, given a vector anysize
%
% input: 
% indx = vector of positions
% n = matrix to reshape (only needed for size)!!!
%
indx = varargin{1};
n = size(varargin{2},1);

% make the vector large as it should
x = indx;

v = 1:(n*(n-1)/2);

[a,b] = find(v == x);
ve = zeros(1,(n*(n-1)/2));
ve(b) = a;

% do the matrix
A = triu(ones(n),1);
A(A ~=0 ) = ve;

% localize the two nodes
[a1, b1, val] = find(A);

