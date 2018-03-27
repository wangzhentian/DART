function o = convertCor(xi,yi,Nx,Ny,pixelSize,theta)
%双线性差值的坐标变换
%默认pixelSize == 1
%
o = zeros(4,3);
%第一个点
%矩阵坐标到像素空间坐标的变换
x = (-Nx/2 + xi - 0.5);
y = (Ny/2 - yi + 0.5);
%变换前的空间坐标
%正旋
t = cos(theta)*x-sin(theta)*y;
s = sin(theta)*x+cos(theta)*y;
%反旋
% t = cos(theta)*x+sin(theta)*y;
% s = -sin(theta)*x+cos(theta)*y;

%空间坐标到矩阵坐标的变换
mx = Nx/2+0.5+t; 
my = Ny/2+0.5-s;

mxI = floor(mx);
mxR = mx-mxI;
myI = floor(my);
myR = my-myI;

% 可由原图像中坐标为 (i,j)、(i+1,j)、(i,j+1)、(i+1,j+1)所对应的周围四个像素的值决定，即：
%f(i+u,j+v) = (1-u)(1-v)f(i,j) + (1-u)vf(i,j+1) + u(1-v)f(i+1,j) + uvf(i+1,j+1)
if mxI>0 && mxI+1 <= Nx && myI > 0 && myI+1 <= Ny
    o(1,:) = [ myI   mxI   (1-mxR)*(1-myR)];
    o(2,:) = [ myI+1 mxI+1 mxR*myR];
    o(3,:) = [ myI+1 mxI   (1-mxR)*myR];
    o(4,:) = [ myI   mxI+1 mxR*(1-myR)];
    %n2(yi,xi) = (1-mxR)*(1-myR)*n(myI,mxI) + mxR*myR*n(myI+1,mxI+1) + (1-mxR)*myR*n(myI+1,mxI) + mxR*(1-myR)*n(myI,mxI+1);
end
