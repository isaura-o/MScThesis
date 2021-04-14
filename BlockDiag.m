function A = BlockDiag(varargin)
% inputs:
% varargin{1} = networks input (18 x 18 with Ang nodes in the
% positions 9 and 10)
% 
% 
% outputs:
% block diagonal network (16x16)
%
%
%
mat = varargin{1};

blockD = mat(1:8,1:8);
blockF = mat(11:18,11:18);
block = blkdiag(blockD, blockF);

A = block;
    
end

