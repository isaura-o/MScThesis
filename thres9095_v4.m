function [val99] = thres9095_v4(Xs,th,dim)
% 
% Threshold matrix using certain value as threshold (th) levels through the surrogates
% adjacency matrices. Works for CC and PC. GLOBAL THRESHOLD
%
%
% input:
% Xs = CC or PC adjacency matrix from surrogates series (ROI x ROI x 1000 x Subject)
% th = threshold in values 0 to 1 (not percent!)
% dim =
%       1 => concatenated networks
%       2 => subject networks
%       3 => subject recording
% 
% output:
% val99 = matrix used for threshold matrices.
% 


if (dim == 1)

    tmp1 = Xs;
    
elseif (dim == 2)
    
    tmp1 = reshape(Xs, size(Xs,1),size(Xs,2), []);
    
elseif (dim == 3)

    tmp1 = reshape(Xs,size(Xs,1), size(Xs,2),[]);
elseif dim == 4
    tmp1 = reshape(Xs,size(Xs,1),size(Xs,2),[]);
    
else
    fprintf('Write correct dimension: 1 or 2')
end

val99up = zeros(size(Xs,1));
for i = 1:size(Xs,1)
    for j = i+1:size(Xs,1)
       newsur = tmp1(i,j,:);
       [va, vsr]  = ecdf(newsur(:));
       val99up(i,j) = vsr(floor(th*numel(vsr)));
    end
end

val99 = val99up'+val99up;

end