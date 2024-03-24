clear;
clc;
a=1;
dx=0.01;
x=0:dx:1;
dt=0.00001;
%%更改t的终止时间
t=0:dt:0.1;
u=zeros(length(x),length(t));
u(:,1)=sin(pi*x);   %初始温度分布为sin(x)
m1=0+0*sin(t);    %左边界条件
m2=0+0*cos(t);    %右边界条件
A=-2*eye(length(x))+diag(ones(1,length(x)-1),1)+diag(ones(1,length(x)-1),-1);   %生成矩阵-1 2 -1的矩阵
for n=1:length(t)-1
    u(:,n+1)=u(:,n)+a^2*dt/dx^2*A*u(:,n);
    u(1,n+1)=m1(n+1);
    u(end,n+1)=m2(n+1);
end
ua=exp(-pi*pi*dt*n)*sin(pi*x);
plot(x,u(:,n),'r',x,ua,'b','LineWidth',1.2);
legend("数值解","解析解");xlabel("x/m");ylabel("Temperature/K");axis([0 1 0 1]);grid on;title("t=0.1s时温度分布");