clear
clc
global A1 A2 A3 ENDA Z HMO DX X KK AK N DH 
N=130;
X0=-4;
XE=1.5;
W=1e5;
E1=2.2e11;
EDA0=0.05;
R=0.05;
US=1.5;
Z=0.68;
Po=1.96e8;
PAI1=0.318309886;
%%求最小油膜厚度
W1=W/(E1*R);
PH=E1*sqrt(0.5*W1/pi);
A1=log(EDA0)+9.67;
A2=PH/Po;
A3=0.59/(PH*1e-9);
B=4*PH*R/E1;
ALFA=Z*A1/Po;
G=ALFA*E1;
U=EDA0*US/(2*E1*R);
CC1=sqrt(2*U);
AM=2*pi*(PH/E1)^2/CC1;
ENDA=3*(pi/AM)^2/8;
HMO=1.6*(R/B)^2*G^0.6*U^0.7*W1^(-0.13);
%% 弹性变形参数
AK=[];
for i=0:N
    AK(i+1)=(i+0.5)*(log(abs(i+0.5))-1)-(i-0.5)*(log(abs(i-0.5))-1);
end
MK=1;
DX=(XE-X0)/(N-1);
X=[];P=[];
for k=1:N
    X(k)=X0+(k-1)*DX;
    if abs(X(k))>=1
        P(k)=0;
    else 
        P(k)=sqrt(1-X(k)^2);
    end
end
%%计算膜厚、密度和粘度弹性变形
KK=0;
[output1,HOO]= HREE( P ,0);
H=output1(1,:);
EDA=output1(2,:);
RO=output1(3,:);
EPS=output1(4,:);
ht=HOO;
%%计算前一次压力迭代值
POLD=[];
for kk=1:N
    POLD(kk)=P(kk);
end
ERP=1;MK=1;
while ERP>1e-5 && DH>1e-6
   KK=19;
   for k=1:KK
       D2=0.5*(EPS(1)+EPS(2));
       D3=0.5*(EPS(2)+EPS(3));
       for ii=2:N-1
           D1=D2;
           D2=D3;
           if ii~=N-1
               D3=0.5*(EPS(ii+1)+EPS(ii+2));
           end
           D8=RO(ii)*AK(1)*PAI1;
           D9=RO(ii-1)*AK(2)*PAI1;
           D10=1/(D1+D2+(D9-D8)*DX);
           D11=D1*P(ii-1)+D2*P(ii+1);
           D12=(RO(ii)*H(ii)-RO(ii-1)*H(ii-1)+(D8-D9)*P(ii))*DX;
           P(ii)=(D11-D12)*D10;
           if P(ii)<0
               P(ii)=0;
           end
           [output2,HOO]= HREE( P ,ht);
           H=output2(1,:);
           EDA=output2(2,:);
           RO=output2(3,:);
           EPS=output2(4,:);
           ht=HOO;
       end
       MK=MK+1;
       SD=0;
       SUM=0;
       for kk=1:N
          SD=SD+abs(P(kk)-POLD(kk)); 
          POLD(kk)=P(kk);
          SUM=SUM+P(kk); 
          ERP=SD/SUM;
       end
       if MK>=50
           MK=1;
           DH=0.5*DH;
       end
   end
end
   if DH<=1e-6
       disp('压力不收敛')
   end
   H2=1e3;
   P2=0;
   for jj=1:N
       if H(jj)<H2
           H2=H(jj);
       end
       if P(jj)<P2
           P2=P(jj);
       end
   end
   H3=H2*B^2/R;
   P3=P2*PH;
   plot(X,H);
   hold on
   plot(X,P);
   axis([-4 1.5 0 1.5])

           








