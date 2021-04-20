 function [a, prec, rec] = Variability(varargin) 
%
% input:
%   1. Concatenate network => arg{1}
%   2. Subject networks => arg{2}
%   3. Degree of the concatenated network
%   4. title for the graph
%
%
%colormy = load('Documents/Isaura/colorsForMatx.dat');
%
A = varargin{1};
Y = varargin{2};
x = varargin{3};
tit = varargin{4}; %for graph

prec = zeros(1,size(Y,3)); rec = zeros(1,size(Y,3));
for i = 1:size(Y,3)
    [ prec(i), rec(i), ~, ~ ] = assessPerformance( A, Y(:,:,i) );
end


% For figure:
  cen = {(1/x):(1/x):1,(1/x):(1/x):1};
  [a, b] = hist3([prec', rec'], cen);
%  figure; 
%  pcolor(a/numel(rec)); title(tit); xlabel('recall'); ylabel('precision'); 
%  set(gca,'xtick', 0.1:(1/x):1); set(gca, 'ytick', 0.1:(1/x):1); 
%  set(gca,'yticklabel',b{1}); set(gca,'xticklabel',b{2});
%  
%  colorbar('location','Manual', 'position', [0.01 0.1 0.02 0.81]);
%  colormap(colormy)


end