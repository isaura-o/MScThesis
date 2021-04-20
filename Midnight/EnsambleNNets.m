function [netMax, val] = EnsambleNNets(Ypc,L,N,n)
% Calls SmaxCMI to infer the networks using MECMI method. 
% Since Maximum Entropy is a local maximum we permute N times the nodes to
% obtain the maximum Maximum Entropy, which will give the network of the
% system.
% 
%
% Inputs
% Ypc = time series (nodes x TP).
% L = number of nodes.
% N = number of permutations (number of networks in the ensamble. Default at 10)
% n = number of bins (data binnerized in 2, 3 or 4 states).
%
% Outputs
% netMax = network that has the max in maximum entropy of the ensamble.
% val = Value of the maximum entropy.
%
% NOTE: SmaxCMI is the function that infers the network by using MECMI
% algorithm which was written by Elliot Martin, thus is not in the folder.
%
pMat = zeros(L,L,N);
Smax = zeros(1,N);
nodes = zeros(N,L); 
for i = 1:N
    % permute the nodes.
    nodes(i,:) = randperm(L,L);
    Y = Ypc(nodes(i,:),:);
    % infers the network.
    % Calculates Smax and optimizes the eq system with Smax as boundary condition to get the CMI values.
    [Smax(:,i), pMat(:,:,i)] = SmaxCMI(Y,n); 
    
end    

% Sort the nodes in the order that they are on: DMN, FPN, or DMN+FPN.
[~, indx] = sort(nodes,2);
newMat = zeros(L,L,N);
for i = 1:N
    tmp = pMat(:,:,i);
    newMat(:,:,i) = tmp(indx(i,:), indx(i,:));
end

% Maximum value of Smax
[val, ind] = max(Smax);
% Extract the network that has this Smax.
netMax = newMat(:,:,ind);

end