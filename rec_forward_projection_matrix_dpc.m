%In line DART iterative algorithm.
%����Ҫ�ڴ棬���Դ���large scale problem.
%�����м���ϵ������

Nb = 512;
Na = 360;
Nx = Nb;
Orbit = 180;

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

Nb = Nb2;
Nx = Nb;


%few view limit
Dangle = Orbit/Na;
Ra = 0:Dangle:(Orbit-Dangle);

pixelSize = 1;

%% here I have to generate a sparse matrix
% forward projection matrix
%P = zeros(Nb*numel(Ra),Nx*Nx);
% create temperary space
xxx = [];
yyy = [];
vvv = [];
ddd = [];

Ny = Nb*numel(Ra);

parpool('local',16);

%for ai = 0:Dangle:Na-1
    %theta = pi/180*ai;
parfor ai = 1:numel(Ra)
    tic
    theta = pi/180*Ra(ai);
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
            elseif ti == Nb %ͼ�����ұߵ�һ��,��Ϊ���ǵ����ϵ�ͼ��û�������
                %��һ����ľ�������
                x1 = ti-1;
                y1 = yi;
                %�ڶ�����ľ�������
                x2 = ti;
                y2 = yi; 
            else
                %��һ����ľ�������
                x1 = ti-1;
                y1 = yi;
                %�ڶ�����ľ�������
                x2 = ti;
                y2 = yi; 
            end

            o1 = convertCor(x1,y1,Nx,Nx,pixelSize,theta);
            o2 = convertCor(x2,y2,Nx,Nx,pixelSize,theta);

            %o1 ��֣�����
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
        vp = (ai-1) * Nb + ti;
        ids = find(C~=0);
        xxx = [xxx ids];
        yyy = [yyy ones(size(ids))*vp];
        vvv = [vvv C(ids)];

    end%detector
    toc
end%angle

P = sparse(yyy,xxx,vvv,Ny,Nx*Nx,numel(vvv));
delete(gcp('nocreate')); % stop the para pool
