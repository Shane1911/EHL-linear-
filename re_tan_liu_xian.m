clear
clc
N=130;
aa=-4;bb=1.5;
v1=0.26;v2=0.26;
E=2.07e11;%����ģ��
X=linspace(aa,bb,N);%lg1h
DX=(bb-aa)/N;
H0=0;
H=X.^2/2;
T0=303; 
p_0=1.96e8;%ϵ��
yita0=0.03;%p=0 T=T0ʱ�󻬼�ճ��
w=1.76e5;%��λ�����غ�
W=pi/2;
R=0.02;%�����뾶
E=2.11e11;%����ģ��
D=-0.00065;%����ϵ��  
b=sqrt(8*w*R/(pi*E));%�߽Ӵ������
p0=1.96e8;%���ѹ��
ph=sqrt(w*E/(2*pi*R*(1-v1^2)));%���ȽӴ�Ӧ������Բ�ֲ�
T=ones(1,N);%�����ٻ���ʼ�¶�
c=1e-3;
us=1.5;%
lamda=12*yita0*us*R^2/(b^3*ph);
a0=((0.5)*(log(abs(0.5))-1)-(-0.5)*(log(abs(-0.5))-1))/DX-log(DX);%�ڵ�0����ϵ��
a1=((1+0.5)*(log(abs(1+0.5))-1)-(1-0.5)*(log(abs(1-0.5))-1))/DX-log(DX);%�ڵ�1����ϵ��
% Hmin=1.6*(R/b)*(0.49*E)^0.6*(us*yita0/E/R)^0.7/((w/E/R)^0.13);%�����ٻ�����СĤ��
hmin=1.6*(2.2e-8)*(us*yita0)^0.7*R^0.43*E^0.03/(w^0.13)
Hmin=hmin*R/b^2;
detaH0=0.005*Hmin;%�𲽼�С����
PP=zeros(1,N);
k=1;
kk=1;
for i=1:N;
    if abs(X(i))>1;
        P(i)=0;
    else
        P(i)=sqrt(1-X(i)^2);
    end
end %P:��ʼѹ��


%%
zaihe=sum(P)*DX;
detaw=pi/2-zaihe;



%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
while max(abs((PP-P)./P))>c;  
%% Ĥ�� 

vi=zeros(1,N);%���Ա���
 for i=1:N; 
      for j=1:N; 
        I=abs(i-j);  
       Kij(i,j)=(I+0.5)*(log(abs(I+0.5))-1)-(I-0.5)*(log(abs(I-0.5))-1);%���Ա���ϵ����
       vi(i)=vi(i)+(Kij(i,j)+log(DX))*P(j)*DX/pi;
      end
 end
   H=H-vi; 
%% ����ճ�� �ܶ�

YITA=exp((log(yita0)+9.67).*((1+ph.*P./p_0).^(0.68)-1));%niandu p
RO=1+0.6.*P./(1+1.7.*P);
%% ��ŵ���� ������ѹ���ֲ�
PP=P;%������ѹ��  
while 1
 for i=2:N-1; 
ei=RO(i)*H(i)^3/(YITA(i)*12*yita0*us*R^2/(b^3*ph));
ei1=RO(i+1)*H(i+1)^3/(YITA(i+1)*12*yita0*us*R^2/(b^3*ph));
ei_1=RO(i-1)*H(i-1)^3/(YITA(i-1)*12*yita0*us*R^2/(b^3*ph));
ei1_2=(ei+ei1)/2;
ei_1_2=(ei+ei_1)/2;
P(i)=(ei_1_2*P(i-1)+ei1_2*PP(i+1)-DX*PP(i)*(RO(i)*H(i)-RO(i-1)*H(i-1)+(RO(i)*a0/pi-RO(i-1)*a1/pi)))/(ei_1_2+ei1_2+(RO(i-1)*a1/pi-RO(i)*a0/pi)*DX);
 end
%�غ�ƽ�ⷽ��
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
%% �����������¶�
% for i=1:N;
% T(i)
 

end 
plot(X,P,'r');
hold on
plot(X,H,'b');

    
    

    
    



