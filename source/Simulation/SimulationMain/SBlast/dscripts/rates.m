clear;
load error.dat;

[n,m]=size(error);

x=error(:,1);
d=error(:,2);
p=error(:,3);
u=error(:,4);

xx=log10(x);
bd=log10(d);
bp=log10(p);
bu=log10(u);

A=ones(n,2);
A(:,2)=xx;

dbeta=lsqr(A,bd,1);
pbeta=lsqr(A,bp,1);
ubeta=lsqr(A,bu,1);

yd=10.^(A*dbeta);
yp=10.^(A*pbeta);
yu=10.^(A*ubeta);

y=[x, yd, yp, yu];

save -ascii ratefit y;
