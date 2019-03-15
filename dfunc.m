% function f = func( x0,y0,z0,l,m,n,dx,dy,dz,r,d )
function [df,f] = dfunc( dx,dy,dz,d,r )
%FUNC 此处显示有关此函数的摘要
%   此处显示详细说明
% syms x0 y0 z0 l m n 
syms x0 y0 z0 theta phi 
% f = (x0+dx-l*d)^2 + (y0+dy-m*d)^2 +(z0+dz-n*d)^2 -r^2;
f = (x0+dx-sin(theta)*cos(phi)*d)^2 + (y0+dy-sin(theta)*sin(phi)*d)^2 +(z0+dz-cos(theta)*d)^2 -r^2;
% df = [diff(f,'x0') ; diff(f,'y0');diff(f,'z0');...
%     diff(f,'l') ; diff(f,'m');diff(f,'n')];
df = [diff(f,'x0') , diff(f,'y0'),diff(f,'z0'),...
    diff(f,'theta') , diff(f,'phi')];
end

