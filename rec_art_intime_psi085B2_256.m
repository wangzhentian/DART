%In line DART iterative algorithm.
%����Ҫ�ڴ棬���Դ���large scale problem.
%�����м���ϵ������
clear all;
load sino_psi085B2_256_real180;

Nb = 256;
Na = 180;
Nx = Nb;

%����ԭʼ��ͼ���С
Nbo = Nb;
%������Ҫб�ߵ�����ֵ
Nb2 = Nb;
% Nb2 = ceil(sqrt(2)*Nb);
if mod(Nb2,2) == 1
    %����
    Nb2 = Nb2 +1;
else
    %ż��
end
%��Ҫ��չ����ͼ
%����ͼ�з����ǽǶ�,�з���������
sino2 = zeros(Na,Nb2);
t = (Nb2-Nb)/2;
sino2(:,t+1:t+Nb) = sino;
sino = sino2;
Nb = Nb2;
Nx = Nb;

%TV parameters
isTV  = 0;
nIter = 20;
e   = 1e-8;
dt    = 1e-5;


%ART ��������
Niter = 200;

% sino2 = zeros(Na,Nb);
% sino2(:,(Nb/2-31):(Nb/2-32+64)) = sino;
% y = sino2';

y = sino';

%few view limit
Dangle = 10;
Ra = 1:Dangle:180;
y = y(:,Ra);

y = -y(:);

pixelSize = 1;
%image,initial value
x = ones(Nx*Nx,1)*1e-6;
xs = zeros(size(x,1)*size(x,2),Niter);

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
                Dx = 1;% ��pixelSize�й�
                if ti == 1 %ͼ������ߵ�һ��,��Ϊ���ǵ����ϵ�ͼ��û�������
                    %��һ����ľ�������
                    x1 = ti;
                    y1 = yi;
                    %�ڶ�����ľ�������
                    x2 = ti+1;
                    y2 = yi;   
                    Dx = 1;
                elseif ti == Nb %ͼ�����ұߵ�һ��,��Ϊ���ǵ����ϵ�ͼ��û�������
                    %��һ����ľ�������
                    x1 = ti-1;
                    y1 = yi;
                    %�ڶ�����ľ�������
                    x2 = ti;
                    y2 = yi; 
                    Dx = 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 elseif ti == Nb/2 %ͼ���м������,��Ϊ�м���������㣬����һ��
%                     %��һ����ľ�������
%                     x1 = ti-1;
%                     y1 = yi;
%                     %�ڶ�����ľ�������
%                     x2 = ti+1;
%                     y2 = yi; 
%                     Dx = 1;    
%                 elseif ti == Nb/2+1 %ͼ���м������,��Ϊ�м���������㣬����һ��
%                     %��һ����ľ�������
%                     x1 = ti-1;
%                     y1 = yi;
%                     %�ڶ�����ľ�������
%                     x2 = ti+1;
%                     y2 = yi; 
%                     Dx = 1;        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% add for test 
%                 elseif ti <= Nb/2 %���ĵ�����
%                     %��һ����ľ�������
%                     x1 = ti-1;
%                     y1 = yi;
%                     %�ڶ�����ľ�������
%                     x2 = ti;
%                     y2 = yi; 
%                     Dx = 1;              
%                 elseif ti > Nb/2 %���ĵ�����
%                     %��һ����ľ�������
%                     x1 = ti;
%                     y1 = yi;
%                     %�ڶ�����ľ�������
%                     x2 = ti+1;
%                     y2 = yi; 
%                     Dx = 1;     
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                else
                    %��һ����ľ�������
                    x1 = ti-1;
                    y1 = yi;
                    %�ڶ�����ľ�������
                    x2 = ti;
                    y2 = yi; 
                    Dx = 1;
                end

                o1 = convertCor(x1,y1,Nx,Nx,pixelSize,theta);
                o2 = convertCor(x2,y2,Nx,Nx,pixelSize,theta);

                %o2 ��֣�����
                for in = 1:4
                    yy = o1(in,1);
                    xx = o1(in,2);
                    w  = o1(in,3);
                    if yy ~= 0 && xx ~=0 && yy <= Nx && xx <= Nx
                        %��ά����ת��Ϊ��������
                        vx = (xx-1) * Nx + yy;
                        C(1,vx) = C(1, vx) + w/Dx;
                    end
                end
                %o2 ��֣�����
                for in = 1:4
                    yy = o2(in,1);
                    xx = o2(in,2);
                    w  = o2(in,3);
                    if yy ~= 0 && xx ~=0 && yy <= Nx && xx <= Nx
                        %��ά����ת��Ϊ��������
                        vx = (xx-1) * Nx + yy;
                        C(1,vx) = C(1, vx) - w/Dx;
                    end
                end  
                
            end%ray  
            
            %�ڼ�������
            vp = ai/Dangle * Nb + ti;
            
            %һ�����ߵĵ�������
            x0 = x;
            a = C;
            yi = a*x;
            c = a*a';
            if c ~= 0
              x = x + a'.*(y(vp)-yi)/c;
            end          
            
            %regurization
            x(x<0) = 0;
            
        end%detector
    end%angle
    
    %TV iteration
    if (isTV && iter > 0)
        x = x*1e8;
        x2 = x2*1e8;
        d=norm(x-x2);
        x = reshape(x,Nx,Nx);
        x = tv_pan(x,nIter,dt,d,e);
        x = x(:);
        x = x*1e-8;
        x2 = x2*1e-8;
    end
    cov(iter) = sum(sum(abs(x-x2)));
    
    xs(:,iter) = x;
    toc
end %iter

x = reshape(x,Nx,Nx);
x = x';
figure;imshow(x,[]);

% x = x(t+1:t+Nbo,t+1:t+Nbo);
% figure;imshow(x,[]);

%phantom
% N = 64;
% fid = fopen('exact\0.dat','r');
% img = fread(fid,[N,N],'double');
% fclose(fid);
% img = img';
% figure;imshow(img,[]);
% figure;plot(img(36,:));hold on;plot(x(36,:),'r');