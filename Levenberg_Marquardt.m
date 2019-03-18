% Levenberg_Marquardt
format long
data = load('sixpoint12152119.txt');
% [l m n ] = [1 1 1]/sqrt(3)
% [x0 y0 z0] = [30 20 10];
% load('point_randd.mat');
% data = daxyz;

% Levenberg-Marquardt
x = data(:,1);
y = data(:,2);
z = data(:,3);
d = data(:,4);
R = 50.7698/2;

x_0 = [0,0,0,0,-pi];%������ֵ
% #######################
% ����x0 lb&ub 
% #######################
dx = data(:,1)-data(1,1);
dy = data(:,2)-data(1,2);
dz = data(:,3)-data(1,3);
dxyz = [dx,dy,dz];
DfT=sym(zeros(length(x_0),length(x)));
F=DfT(:,1);
% �������ӳ�ʼ�� �Լ��������Ӹı���
mu = 20;
v = 1.5;
IDFF = eye(5);
%
for i=1:length(x)
    [df,f] = dfunc( dxyz(i,1),dxyz(i,2),dxyz(i,3),data(i,4),R );
    %���ż����� A��������Hermitת�� A.'������ת��
    DfT(:,i)=df.';% ��������DfT��ʽһ�� ��һ������������С����������·�����
    F(i) = f;%Ϊ������
end
xk = x_0';
digits(8);
% �ڶ�������ѭ�� ��һ������
while(1)
    xk
    DF = subs(DfT,{'x0','y0','z0','theta','phi'},{xk(1),xk(2),xk(3),xk(4),xk(5)});
    FF = subs(F,{'x0','y0','z0','theta','phi'},{xk(1),xk(2),xk(3),xk(4),xk(5)});
    FFF = double(vpa(FF));
    DFF = double(vpa(DF));
    Gx = (inv(DFF*DFF'+mu*IDFF))*DFF*FFF;
    xk1 = xk-Gx;
    % �˴������ܶ��
    FF1 = subs(F,{'x0','y0','z0','theta','phi'},{xk1(1),xk1(2),xk1(3),xk1(4),xk1(5)});
    FFF1 = double(vpa(FF1));
    if(abs(FFF)<=abs(FFF1))
        mu = mu*v;
    else
        xk = xk1;
        mu = mu/v;
%         if(abs(DFF*FFF)<=1e-2)
        if(abs(Gx)<=1e-2)
            break;
        end
    end
end

xk1
theta=xk1(4);phi=xk1(5);
l=sin(theta)*cos(phi)
m=sin(theta)*sin(phi)
n=cos(theta)