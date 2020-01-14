function [ Y h] = HREE( P,ht)
global X A1 A2 A3 ENDA Z HMO DX KK N DH 
Go=pi/2;
if KK==0
    HOO=0;
else 
    HOO=ht;
end
W1=0;
for i=1:N
   W1=W1+P(i);
end
c3=(DX*W1)/Go;
DW=1-c3;
V=VI(P,N);
HMIN=1e3;
H=[];
for k=1:N
    HO=0.5*X(k)^2+V(k);
    if HO<HMIN
        HMIN=HO;
    end
    H(k)=HO;
end
if KK==0
  KK=1;
  DH=0.005*HMO;
  HOO=-HMIN+HMO;
end
if DW<0
    HOO=HOO+DH;
else 
    HOO=HOO-DH;
end
EDA=[];RO=[];EPS=[];
for kk=1:N
    H(kk)=H(kk)+HOO;
    EDA(kk)=exp(A1*(-1+(1+A2*P(kk))^Z));
    RO(kk)=(A3+1.35*P(kk))/(A3+P(kk));
    EPS(kk)=RO(kk)*H(kk)^3/(ENDA*EDA(kk));
end
 Y(1,:)=H;
 Y(2,:)=EDA;
 Y(3,:)= RO;
 Y(4,:)= EPS;
 h=HOO;
end

