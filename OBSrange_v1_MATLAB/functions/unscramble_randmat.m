function [ matout ] = unscramble_randmat( mat, indxs )
% [ matout ] = unscramble_randmat( mat, indxs )
%  
% Takes a matrix of values with indexes 1:Npara which have been randomized and
% unscrambles them. 
% 
% IN:
% mat : matrix with random indexes
% indxs : matrix containing random indexes contains the indexed values of each
%
[Nind ,Mperms] = size(indxs);
matout = zeros(size(indxs));
for iind = 1:Nind
    matout(iind,:) = mat(indxs==iind);
end


end

