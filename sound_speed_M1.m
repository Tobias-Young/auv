%仰涛，2017/11/14.kraken 模型。
%%simulation of the ***
%%2017-11-15 修改了**部分

clc;clear all;close all;
%%定义变量，可以定义类型、精度，
H=3000;   R=5000;   %水深3000米，宽度5000米
C0=1500; pho=1025; %基准声速和水密度
freq=15;  %频率
N=300; N1=500; %深度分层和水平分层层数
h=H/N; r=R/N1; %分层水柱厚度
SD=3000; sh=SD/h; %声源深度和所在层数
%%定义边界、初始条件等，简单函数
z=((1:N) -1/2)*h;    %每一层深度，对声速数据进行插值，分N=300层，每层深度1米。
x=((1:N1) -1/2)*r; 
c=C0*[1+0.00737*(2*(z-1300)/1300-1+exp(-2*(z-1300)/1300))]; %c 是声速，盐度、温度、密度、浑浊度等影响因素
%cc=interp1((1:H),c,z,'spline');  %三次样条插值,x,y为分300层后的声速剖面数据
%cc=c(z);  %取出中间300层的声速剖面数据
figure(1); 
plot(c,z,'k');  %画声速剖面图。
set(gca,'YDir','reverse');
title('声速剖面图');
legend('插值前声速','插值后声速')
xlabel('c（m/s）'); ylabel('d(m)');
%%矩阵A的求取
w=2*pi *freq;  %声源角频率
for i=1:N
    A(i,i)= -2+(h*w)^2/c(i)^2;
end
for i=1:N-2
    A(i,i+1)=1;
    A(i+1,i)=1;
end
A(N-1,N)=1; 
A(N,N-1)=2;
%%
%求取A的特征向量V和特征值D
[V,D]=eig(A);
for i=1:N
    k(i)=(D(i,i)/(h^2))^0.5;   %k(i)为第i个波数
end
%%
% 求声压p
for j=1:N1    %j为水平距离
for i=1:N  %i为层数
pp=0;
for ii=1: N     %ii为模函数的数量
     pp=pp+(V(sh,ii)*V(i,ii))*exp(sqrt(-1)*k(ii)*j*r)/(k(ii)^0.5); %声源深度V（？,ii）
end
p(j,i)=(sqrt(-1)/(pho*((8*pi*j*r)^0.5)))*(exp(-sqrt(-1)*pi/4))*(pp);
end
end
k0=w/c(sh);  %c(SD)，声源处的波数
p0=exp(sqrt(-1)*1*k0)/(4*pi*1);  %自由声场中1米处的声压
spl=20*log10(abs(p)/0.00002);    %spl为声压级
TL=-20*log10(abs(p./p0));   %传播损失
save('kraken-m1','H','R','N','N1','SD','sh','freq','pho','c','z','x','TL','p','spl');
pause
%%
figure;
mesh(r,z,(TL)');
view(0,90);
set(gca,'YDir','reverse');
title('传播损失图');
xlabel('r/m');
ylabel('z/m');
colorbar;
%%
%求某一接收深度处的传播损失
figure;
TLL=TL(r,300);
plot(r,TLL);
title('3000米深处传播损失图');
xlabel('r/m');
ylabel('z/m');
%%
%绘制声压图
figure(3);
mesh(x,z,(spl)');
view(0,90);
set(gca,'YDir','reverse');
title('声压图');
xlabel('z/m');
ylabel('x/m');
c=colorbar;
c.Label.String = '声压级/dB';
%绘制声压水平剖面图
figure(4);
plot(x,spl(x,300))
xlabel('r/m');
ylabel('声压级/dB');
title('声压水平剖面(z=2995m)');
%绘制声压垂直剖面图
figure(5);
plot(z,(spl(1,:))')
xlabel('z/m');
ylabel('声压级/dB');
title('声压垂直剖面(r=0m)');
