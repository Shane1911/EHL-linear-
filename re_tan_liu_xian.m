clear
clc
N=130;
aa=-4;bb=1.5;
v1=0.26;v2=0.26;
E=2.07e11;%弹性模量
X=linspace(aa,bb,N);%lg1h
DX=(bb-aa)/N;
H0=0;
H=X.^2/2;
T0=303; 
p_0=1.96e8;%系数
yita0=0.03;%p=0 T=T0时润滑剂粘度
w=1.76e5;%单位长度载荷
W=pi/2;
R=0.02;%当量半径
E=2.11e11;%弹性模量
D=-0.00065;%密温系数  
b=sqrt(8*w*R/(pi*E));%线接触区半宽
p0=1.96e8;%入口压力
ph=sqrt(w*E/(2*pi*R*(1-v1^2)));%赫兹接触应力按椭圆分布
T=ones(1,N);%无量纲化初始温度
c=1e-3;
us=1.5;%
lamda=12*yita0*us*R^2/(b^3*ph);
a0=((0.5)*(log(abs(0.5))-1)-(-0.5)*(log(abs(-0.5))-1))/DX-log(DX);%节点0弹性系数
a1=((1+0.5)*(log(abs(1+0.5))-1)-(1-0.5)*(log(abs(1-0.5))-1))/DX-log(DX);%节点1弹性系数
% Hmin=1.6*(R/b)*(0.49*E)^0.6*(us*yita0/E/R)^0.7/((w/E/R)^0.13);%无量纲化的最小膜厚
hmin=1.6*(2.2e-8)*(us*yita0)^0.7*R^0.43*E^0.03/(w^0.13)
Hmin=hmin*R/b^2;
detaH0=0.005*Hmin;%逐步减小步长
PP=zeros(1,N);
k=1;
kk=1;
for i=1:N;
    if abs(X(i))>1;
        P(i)=0;
    else
        P(i)=sqrt(1-X(i)^2);
    end
end %P:初始压力


%%
zaihe=sum(P)*DX;
detaw=pi/2-zaihe;



%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
while max(abs((PP-P)./P))>c;  
%% 膜厚 

vi=zeros(1,N);%弹性变形
 for i=1:N; 
      for j=1:N; 
        I=abs(i-j);  
       Kij(i,j)=(I+0.5)*(log(abs(I+0.5))-1)-(I-0.5)*(log(abs(I-0.5))-1);%弹性变形系数项
       vi(i)=vi(i)+(Kij(i,j)+log(DX))*P(j)*DX/pi;
      end
 end
   H=H-vi; 
%% 计算粘度 密度

YITA=exp((log(yita0)+9.67).*((1+ph.*P./p_0).^(0.68)-1));%niandu p
RO=1+0.6.*P./(1+1.7.*P);
%% 雷诺方程 计算新压力分布
PP=P;%保留老压力  
while 1
 for i=2:N-1; 
ei=RO(i)*H(i)^3/(YITA(i)*12*yita0*us*R^2/(b^3*ph));
ei1=RO(i+1)*H(i+1)^3/(YITA(i+1)*12*yita0*us*R^2/(b^3*ph));
ei_1=RO(i-1)*H(i-1)^3/(YITA(i-1)*12*yita0*us*R^2/(b^3*ph));
ei1_2=(ei+ei1)/2;
ei_1_2=(ei+ei_1)/2;
P(i)=(ei_1_2*P(i-1)+ei1_2*PP(i+1)-DX*PP(i)*(RO(i)*H(i)-RO(i-1)*H(i-1)+(RO(i)*a0/pi-RO(i-1)*a1/pi)))/(ei_1_2+ei1_2+(RO(i-1)*a1/pi-RO(i)*a0/pi)*DX);
 end
%载荷平衡方程
detaw=pi/2-sum(P)*DX;
if abs(detaw)<1e-6
    break
end
     if detaw>0
        H=H+detaH0;
       else if detaw<0
        H=H-detaH0;
            end
     end  
kk=kk+1;
end
k=k+1;
%% 能量方程算温度
% for i=1:N;
% T(i)
 

end 
plot(X,P,'r');
hold on
plot(X,H,'b');

    
    

    
    



