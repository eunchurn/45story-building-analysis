clear all; clc; close all;

n=3;
m1=100;  % kg
k1=500000; % N/m
c1=100000; 
K=zeros(n+1,n+1);
C=zeros(n+1,n+1);
for kk=1:n
    M(kk,kk)=m1;
    K(kk:kk+1,kk:kk+1)=K(kk:kk+1,kk:kk+1)+[k1 -k1;-k1 k1];
    C(kk:kk+1,kk:kk+1)=C(kk:kk+1,kk:kk+1)+[c1 -c1;-c1 c1];
end
K=K(2:n+1,2:n+1);
C=C(2:n+1,2:n+1);

[phi,lambda]=eig(K,M);
fn=diag(lambda)/(2*pi);
As=[zeros(n,n) eye(n,n); -inv(M)*K -inv(M)*C];
Bs=[zeros(n,1);-ones(n,1)];
Fs=1000;
in=load('elcacc1.txt')/1000*9.8;
[Ad,Bd]=c2d(As,Bs,1/Fs);

Xo=zeros(2*n,length(in));
for kk=1:length(in)
    Xo(:,kk+1)=Ad*Xo(:,kk)+Bd*in(kk);
end
Cs=As(n+1:2*n,:);
y=Cs*Xo;

figure,plot(Xo(1:n,:)')
figure,bar(max(abs(diff(Xo(1:n,:))')))
figure,plot(diff(Xo(1:n,:)'))