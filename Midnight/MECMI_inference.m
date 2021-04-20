function [netM, val] = MECMI_inference(varargin)
%
% MECMI inference
%
% 1. Concatenate (A = concatenated data set) 
% 2. Subjects (B = subject data set), (rhos = non concatenated surrogate networks)
% 3. Threshold (needed input of matrix of thresholds: thr)
%
% Inputs:
%varargin{1} is the Time - series.
%varargin{2} is the number of nodes we want to infer.
%varargin{3} is the size of the ensamble or number of permutations we want
%                   to do in the nodes (Default is 10 permutations).
% n = number of bins (data binnerized in 2, 3 or 4 states).

%varargin{5} is the dimension of the data: 
%               set it to 2 for concatenated all
%               set it to 3 for subject.
%               set it to 4 for recording-subject.
% Outputs: 
%       netM => MECMI network.
%       val => Smax value
%
A = varargin{1};
L = varargin{2}; 
N = varargin{3};
n = varargin{4};

if varargin{5} == 2
    % concatenated
    Ypc = A - mean(A,2); 
    [netM, val] = EnsambleNNets(Ypc,L,N,n);
    
    
elseif varargin{5} == 3
    % subject (S) network
    A = permute(A, [1 4 2 3]); 
    Ac = reshape(A, size(A,1), size(A,2), []);
    Ypc = Ac - mean(Ac, 3);
    Ypc = permute(Ypc, [1 3 2]);
    
    netM = zeros(size(Ypc,1), size(Ypc,1), size(Ypc,3));
    val = zeros(N,size(Ypc,3));
    for i = 1:size(Ypc,3)
        tmp = Ypc(:,:,i);
        [netM(:,:,i), val(:,i)] = EnsambleNNets(tmp,L,N,n);
    end
       
elseif varargin{5} == 4
    % recording-subject (RS) network
    Ypc = A;
    
    netM = zeros(size(Ypc,1), size(Ypc,1), size(Ypc,3),size(Ypc,4));
    val = zeros(N,size(Ypc,3),size(Ypc,4));
    for i = 1:size(Ypc,3)
        for j = 1:size(Ypc,4)
            tmp = Ypc(:,:,i,j);
            [netM(:,:,i,j), val(:,i,j)] = EnsambleNNets(tmp,L,N,n);
        end
    end
       
end


end
