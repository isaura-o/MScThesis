% driving

A = M3ND1 + M3ND2;
B = PCND1 + PCND2;
A1 = M3ND1;
B1 = PCND1;

%A = abs(varargin{1}); % concatenate 1 (2dim mat) x axis
%B = abs(varargin{2}); % concatenate 2 (2dim mat) y axis
%A1 = abs(varargin{3}); % NB network for method 1
%B1 = abs(varargin{4}); % NB network for method 2


% concatenate
tmpC = A(~tril(ones(size(A))));
% most probable network
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

% Now we want to know the value and the position in the matrix from the links in 
% the vector Y(:,2) that are the values: 22,25,26.
Ycro = (Y(:,1) == 0); Ycrop = Y(:,2).*(Ycro); Ycrop = Ycrop(Ycrop > 0);
Ycro2 = Ycrop(Ycrop <= 20);
z = Ycro2

G = graph(B, 'upper');
B = triu(B,1);
Tr = zeros(1,length(z));
x2 = zeros(1,length(z));
y = zeros(1,length(z));
for l = 1:length(z)
    x = find(pos == z(l)); 
    [a,b] = find(B == val(z(l)))
%    [x2(l), y(l)] = find(B == val(z(l)))


    
    edgxN = neighbors(G,a); % edges of node k
   
   edgsbetween = zeros(size(edgxN,1));
    for i = 1:size(edgxN,1)
        for j = i+1:size(edgxN,1)
            edgsbetween(i,j) = edgecount(G,edgxN(i),edgxN(j)); % number of edges between edges of a node (triangles)
        end
    end
    eneignode = sum(sum(edgsbetween)); %number of triangles on node k
    Tr(l) = eneignode;
end

Tr
%% cluster coefficient
G = graph(B, 'upper')
degs = degree(G)
cl = zeros(size(degs,1),1);
tri = zeros(size(degs,1),1);
for k = 1:size(degs,1)
   edgxN = neighbors(G,k) % edges of node k
   
   edgsbetween = zeros(size(edgxN,1));
    for i = 1:size(edgxN,1)
        for j = i+1:size(edgxN,1)
            edgsbetween(i,j) = edgecount(G,edgxN(i),edgxN(j)); % number of edges between edges of a node (triangles)
        end
    end
    eneignode = sum(sum(edgsbetween)); %number of triangles on node k
   
   
   if degs(k)==1 || degs(k)==0; cl(k)=0; continue; end 
   cl(k) = 2*eneignode/(degs(k)*(degs(k)-1))
   tri(k) = eneignode
end

%% now search if any of these values belong to x2

