function mat=matrixgen(X, Y, coeff, scale)
% generate RBF matrix for collocation X and centers Y
% coeff(1)*(u_xx+u_yy) + coeff(2)*u_x + coeff(3)*u_y + coeff(4)*u
% coeff=[1,0,0,1]


global RBFtype; 
global RBFpar; 
global RBFscale
if nargin<4,
    if max(size(RBFscale))==1,
        scale = zeros(max(size(Y)),1)+RBFscale;
    else
        scale=RBFscale;
    end
end

mat = coeff(1)*laplacekermat(X,Y,scale);

size_coef=length(coeff);
switch size_coef
    case 4%2D coeff(1)*(u_xx+u_yy) + coeff(2)*u_x + coeff(3)*u_y + coeff(4)*u
 
        if coeff(2)~=0 ||  coeff(3)~=0
            tmp = gradkermatX( X,Y,scale);
            mat = mat + coeff(2)*tmp(:,:,1) + coeff(3)*tmp(:,:,2);
            clear tmp
        end
% if coeff(4)~=0
%     if isnan( coeff(4) )
%         global vco4
%         mat = mat + vco4(X)*kermat( X,Y,scale);
%     else
%         mat = mat + coeff(4)*kermat( X,Y,scale);
%     end
% end
    case 5 %3D coeff(1)*(u_xx+u_yy+u_zz) + coeff(2)*u_x + coeff(3)*u_y + coeff(4)*u_z+coeff(5)*u
        if coeff(2)~=0 ||  coeff(3)~=0 ||  coeff(4)~=0
            tmp = gradkermatX( X,Y,scale);
            mat =coeff(4)*tmp(:,:,3);% mat + coeff(2)*tmp(:,:,1) + coeff(3)*tmp(:,:,2)+
            clear tmp
        end
end

if coeff(size_coef)~=0
    if isnan( coeff(size_coef) )
        global vco4
        mat = mat + vco4(X)*kermat( X,Y,scale);
    else
        mat = mat + coeff(size_coef)*kermat(X,Y,scale);
    end
end
%%
function mat=kermat(X, Y, scale)
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for kermat arguments');
end
mat=frbf(DistanceMatrixSquare(X,Y)*diag(scale.^2),0);

function mat=laplacekermat(X, Y, scale)
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for laplacekermat arguments');
end
s=DistanceMatrixSquare(X,Y)*diag(scale.^2);
mat=2*(dx*frbf(s,1)+2*s.*frbf(s,2))*diag(scale.^2);

function mat=gradkermatX(X,Y,scale)
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for gradkermatX arguments');
end
fmat=frbf(DistanceMatrixSquare(X,Y)*diag(scale.^2),1)*diag(scale.^2);
mat=zeros(nx,ny,dx);
for dim=1:dx
    mat(:,:,dim)=2*fmat.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1));
end