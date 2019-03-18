% Levenberg_Marquardt
format long
% data = load('sixpoint12152119.txt');
 load('point_randd.mat');
 data = daxyz;
% Levenberg-Marquardt
x = data(:,1);
y = data(:,2);
z = data(:,3);
d = data(:,4);
R = 50.7698;

x_0 = [0,0,0,0,-pi/4];%迭代初值

dx = data(:,1)-data(1,1);
dy = data(:,2)-data(1,2);
dz = data(:,3)-data(1,3);
dxyz = [dx,dy,dz];
DfT=sym(zeros(length(x_0),length(x)));
F=DfT(:,1);
% 迭代因子初始化 以及迭代因子改变量
mu = 2;
v = 1.5;
IDFF = eye(5);
%
for i=1:length(x)
    [df,f] = dfunc( dxyz(i,1),dxyz(i,2),dxyz(i,3),data(i,4),R );
    %符号计算中 A’求解的是Hermit转置 A.'求解的是转置
    DfT(:,i)=df.';% 与论文上DfT格式一样 《一个求解非线性最小二乘问题的新方法》
    F(i) = f;%为列向量
end
xk = x_0';
digits(8);
% 第二个才能循环 第一组数据
DF = subs(DfT,{'x0','y0','z0','theta','phi'},{xk(1),xk(2),xk(3),xk(4),xk(5)});
FF = subs(F,{'x0','y0','z0','theta','phi'},{xk(1),xk(2),xk(3),xk(4),xk(5)});
FFF = double(vpa(FF));
DFF = double(vpa(DF));
Gx = (inv(DFF*DFF'+mu*IDFF))*DFF*FFF;
xk1 = xk-Gx;
preFFF = FFF; 
if(abs(Gx)<1e-2)
    break;
else
    xk=xk1;
end
% 开始循环
while(1)
    DF = subs(DfT,{'x0','y0','z0','theta','phi'},{xk(1),xk(2),xk(3),xk(4),xk(5)});
    FF = subs(F,{'x0','y0','z0','theta','phi'},{xk(1),xk(2),xk(3),xk(4),xk(5)});
    FFF = double(vpa(FF));
    DFF = double(vpa(DF));
    %     Gx中添加阻尼项u*I
    Gx = (inv(DFF*DFF'+mu*IDFF))*DFF*FFF;
    if(norm(abs(FFF))<=norm(abs(preFFF)))
        xk1 = xk-Gx;
        xk1
    else
        mu = mu*v;
    end

    if(abs(Gx)<1e-3)
        break;
    else
        if(norm(abs(FFF))<=norm(abs(preFFF)))
            preFFF = FFF;
            xk=xk1;
            mu = mu/v;
        end
    end
    
end
xk1
theta=xk1(4);phi=xk1(5);
l=sin(theta)*cos(phi)
m=sin(theta)*sin(phi)
n=cos(theta)