function [T, e1b, e2b] = TrianDrivIndirect(varargin)
% Function for determine "triangles" aka driving/indirect links.
%
% An indirect link can be found applying partial grangers causality. since
% a process A -> B and a process B -> C leads to A -> C. taking into
% account instead that we are using links instead of processes, we can say
% that the nodes in a triangle are the links driven by an indirect link.
%
B = varargin{1};
a = varargin{2};
b = varargin{3};

coorZ = [a,b];
cZUn = zeros(size(coorZ,1), size(coorZ,2));
 p = 1;
   for i=1:(size(coorZ,1)-1)
     r = 1;
       while (((a(r)~= b(i)) || (b(r)~= a(i))) && r < p)
            r = r + 1;  
       end 
        if (r == p) 
         cZUn(r,:) = coorZ(i,:);
         p=p+1;
        end
   end 
    
   cZUn( all(~cZUn,2), : ) = [];
   
   
% take these edges (coordinates) and see how many triangles it has.
% make the matrix logical:
Blog = logical(B);
[e1a, e1b] = find( (Blog(cZUn(:,1),:) == 1)); 
[e2a, e2b] = find((Blog(cZUn(:,2),:) == 1));
T = zeros(1,max(e1a));
for k = 1:max(e1a)
    x1 = find(e1a == k);
    e1 = e1b(x1)';
    x2 = find(e2a == k);
    e2 = e2b(x2)';
    for i =1:length(e1)
        for j = 1:length(e2)
            if (e1(i) == e2(j) )
                T(:,k) = T(:,k) + 1;
            end
        end
    end
end