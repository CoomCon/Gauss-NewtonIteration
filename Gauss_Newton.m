% GaussNewTon
format long
data = load('sixpoint12152119.txt');
% Levenberg-Marquardt
x = data(:,1);
y = data(:,2);
z = data(:,3);
d = data(:,4);
R = 50.7698;

x_0 = [0,0,0,0,-pi];%������ֵ

dx = data(:,1)-data(1,1);
dy = data(:,2)-data(1,2);
dz = data(:,3)-data(1,3);
dxyz = [dx,dy,dz];
DfT=sym(zeros(length(x_0),length(x)));
F=DfT(:,1);
% ��������Է�����
for i=1:length(x)
    [df,f] = dfunc( dxyz(i,1),dxyz(i,2),dxyz(i,3),data(i,4),R );
    %���ż����� A��������Hermitת�� A.'������ת��
    DfT(:,i)=df.';% ��������DfT��ʽһ�� ��һ������������С����������·�����
    F(i) = f;%Ϊ������
end
xk = x_0';
digits(8);
while(1)
    DF = subs(DfT,{'x0','y0','z0','theta','phi'},{xk(1),xk(2),xk(3),xk(4),xk(5)});
    FF = subs(F,{'x0','y0','z0','theta','phi'},{xk(1),xk(2),xk(3),xk(4),xk(5)});
    FFF = double(vpa(FF));
    DFF = double(vpa(DF));
%  *DFF ���� *xk?
    Gx = (inv(DFF*DFF'))*DFF;
    xk1 = xk-Gx*FFF;
    xk1
    if(abs(sum(xk1-xk))<1e-3)
        break;
    else
        xk=xk1;
    end
    
end

