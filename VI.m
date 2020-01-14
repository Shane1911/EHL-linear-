function V=VI(P,N);
global DX AK
V=[];
PAI1=0.318309886;
c=log(DX);
for i=1:N
V(i)=0;
T=V(i);
for j=1:N
    IJ=abs(i-j);
    T=T+(c+AK(IJ+1))*DX*P(j);
end
V(i)=T;
end
for k=1:N
    V(k)=-PAI1*V(k);
end
end

