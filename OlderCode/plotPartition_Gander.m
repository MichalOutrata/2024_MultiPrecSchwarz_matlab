function [ output_args ] = plotPartition_Gander( A, ni, nj )
%PLOTPARTITION Plots a partition of the matrix A
%   Detailed explanation goes here

n = size(A,1);
spy(A);
cni = cumsum(ni);
c0 = {'black','red','black'};
for i=1:length(cni)-1
    col = c0{mod(i-1,3)+1};
    line([0 n],[cni(i)+0.5,cni(i)+0.5],'Color',col);
end
cnj = cumsum(nj);
for j=1:length(cnj)-1
    col = c0{mod(j-1,3)+1};
    line([cnj(j)+0.5,cnj(j)+0.5],[0 n],'Color',col);
end

end

