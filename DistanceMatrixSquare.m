function DM=DistanceMatrixSquare(dsite,ctrs)
[M,s]=size(dsite);[N,s]=size(ctrs);
DM=zeros(M,N);
for d=1:s
    [dr,cc]=ndgrid(dsite(:,d),ctrs(:,d));
     DM=DM+(dr-cc).^2;
end
