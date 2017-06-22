
x=linspace(-1,1,150)';
global RBFtype; 
global RBFscale
global RBFpar; 
RBFtype='g' ;%see frbf.m
RBFpar=1;%
%K=[];
for ep=0.05:0.1:5
    RBFscale= rand(length(x),1)*0+ep;%
    A=matrixgen( x, x.*0,[0,0,0,1], RBFscale);
    %K=[K,A(:,1)];
    figure(1);clf;
    plot(x,A(:,1),'*--')
    title(['\epsilon=' num2str(ep)])
    drawnow();
end
