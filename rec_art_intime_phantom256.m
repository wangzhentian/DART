%In line DART iterative algorithm.
%不需要内存，可以处理large scale problem.
%迭代中计算系数矩阵

load sino_phantom256;

Nb = 256;
Na = 180;
Nx = Nb;

%保存原始的图像大小
Nbo = Nb;
%计算需要斜线的像素值
Nb2 = Nb;
% Nb2 = ceil(sqrt(2)*Nb);
if mod(Nb2,2) == 1
    %奇数
    Nb2 = Nb2 +1;
else
    %偶数
end
%需要扩展正弦图
%正弦图行方向是角度,列方向是像素
sino2 = zeros(Na,Nb2);
t = (Nb2-Nb)/2;
sino2(:,t+1:t+Nb) = sino;
sino = sino2;
Nb = Nb2;
Nx = Nb;

%ART 迭代次数
Niter = 100;
%TV parameters
isTV  = 0;
nIter = 40;
e     = 1e-8;
dt    = 1e-10; 
Dangle = 20;       %few view 角度间隔
isBorder = 1;     %是否加入边界限制条件
BorderWidth = 5;  %限制边界宽度


%保存系数矩阵
%G = sparse(Nb*Na,Nx*Nx);




% sino2 = zeros(Na,Nb);
% sino2(:,(Nb/2-31):(Nb/2-32+64)) = sino;
% y = sino2';

y = sino';


Ra = 1:Dangle:180;
y = y(:,Ra);

y = -y(:);

pixelSize = 1;
%image,initial value
x = ones(Nx*Nx,1)*1e-6;
xs = zeros(size(x,1)*size(x,2),Niter);
x0 = x;

cov = zeros(Niter,1);
for iter = 1:Niter
    tic
    iter
    
     x2 = x;
    for ai = 0:Dangle:Na-1
        theta = pi/180*ai;
        for ti = 1:Nb
            C = zeros(1,Nx*Nx);
            for yi = 1:Nx
                Dx = 1;% 与pixelSize有关
                if ti == 1 %图像最左边的一列,因为考虑到边上的图像没法做差分
                    %第一个点的矩阵坐标
                    x1 = ti;
                    y1 = yi;
                    %第二个点的矩阵坐标
                    x2 = ti+1;
                    y2 = yi;   
                    Dx = 1;
                elseif ti == Nb %图像最右边的一列,因为考虑到边上的图像没法做差分
                    %第一个点的矩阵坐标
                    x1 = ti-1;
                    y1 = yi;
                    %第二个点的矩阵坐标
                    x2 = ti;
                    y2 = yi; 
                    Dx = 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 elseif ti == Nb/2 %图像中间的像素,因为中间出现了亮点，测试一下
%                     %第一个点的矩阵坐标
%                     x1 = ti-1;
%                     y1 = yi;
%                     %第二个点的矩阵坐标
%                     x2 = ti+1;
%                     y2 = yi; 
%                     Dx = 1;    
%                 elseif ti == Nb/2+1 %图像中间的像素,因为中间出现了亮点，测试一下
%                     %第一个点的矩阵坐标
%                     x1 = ti-1;
%                     y1 = yi;
%                     %第二个点的矩阵坐标
%                     x2 = ti+1;
%                     y2 = yi; 
%                     Dx = 1;        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% add for test 
%                 elseif ti <= Nb/2 %中心点的左边
%                     %第一个点的矩阵坐标
%                     x1 = ti-1;
%                     y1 = yi;
%                     %第二个点的矩阵坐标
%                     x2 = ti;
%                     y2 = yi; 
%                     Dx = 1;              
%                 elseif ti > Nb/2 %中心点的左边
%                     %第一个点的矩阵坐标
%                     x1 = ti;
%                     y1 = yi;
%                     %第二个点的矩阵坐标
%                     x2 = ti+1;
%                     y2 = yi; 
%                     Dx = 1;     
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                else
                    %第一个点的矩阵坐标
                    x1 = ti-1;
                    y1 = yi;
                    %第二个点的矩阵坐标
                    x2 = ti;
                    y2 = yi; 
                    Dx = 1;
                end

                o1 = convertCor(x1,y1,Nx,Nx,pixelSize,theta);
                o2 = convertCor(x2,y2,Nx,Nx,pixelSize,theta);

                %o2 差分，正号
                for in = 1:4
                    yy = o1(in,1);
                    xx = o1(in,2);
                    w  = o1(in,3);
                    if yy ~= 0 && xx ~=0 && yy <= Nx && xx <= Nx
                        %二维坐标转换为向量坐标
                        vx = (xx-1) * Nx + yy;
                        C(1,vx) = C(1, vx) + w/Dx;
                    end
                end
                %o2 差分，负号
                for in = 1:4
                    yy = o2(in,1);
                    xx = o2(in,2);
                    w  = o2(in,3);
                    if yy ~= 0 && xx ~=0 && yy <= Nx && xx <= Nx
                        %二维坐标转换为向量坐标
                        vx = (xx-1) * Nx + yy;
                        C(1,vx) = C(1, vx) - w/Dx;
                    end
                end  
                
            end%ray  
            
            %第几条射线
            vp = ai/Dangle * Nb + ti;
            
            %一条射线的迭代过程
            a = C;
            yi = a*x;
            c = a*a';
            if c ~= 0
              x = x + a'.*(y(vp)-yi)/c;
            end          
            
            %regurization
            x(x<0) = 0;

            %加入边界区域限制条件
            if (isBorder)
                x = reshape(x,Nx,Nx);
                x(1:BorderWidth,:) = 0;
                x(end-BorderWidth:end,:) = 0;
                x(:,1:BorderWidth) = 0;
                x(:,end-BorderWidth:end) = 0;
                x = x(:);
            end
                 
        end%detector
    end%angle
    
   
    
    %TV iteration
%     if (isTV && iter > 5)
    if (isTV)    
        x = x*1e8;
        x2 = x2*1e8;
        d=norm(x-x2);
        x = reshape(x,Nx,Nx);
        x = tv_pan(x,nIter,dt,d,e);
        x = x(:);
        x = x*1e-8;
        x2 = x2*1e-8;
    end
    
    xs(:,iter) = x;
    
    cov(iter) = sum(sum(abs(x-x2)));
    %保存一个副本
    x0 = x;
    toc
end %iter

x = reshape(x,Nx,Nx);
x = x';
figure;imshow(x,[]);
