clear all; clc; close all;
% 45층 부산 센텀 리더 건물모델
M=load('m_leader_cen.dat');
C=load('C_leader.dat');
K=load('K_leader.dat');
K2=K;
[r,c]=size(M);
% 자유도
% z = [x1,x2,...,x45,y1,y2,..y45,th1,th2,...,th45]'
% 단위
% ton, meter
n=r;
As=[zeros(n,n) eye(n,n); -inv(M)*K -inv(M)*C];
Bs1=[zeros(n,45);zeros(45,45);inv(M(1:45,1:45));zeros(45,45)];
Bs2=[zeros(n,1);-ones(n*2/3,1);zeros(n*1/3,1)];
Cs=[[eye(n,n) zeros(n,n);zeros(n,n) eye(n,n)];As(n+1:2*n,:)];
Ds=zeros(3*n,1);
C45=Cs(2*n+45,:);
D45=0;
[phi,lambda]=eig(K,M);
Fn=diag(lambda)/(2*pi);
[phix,lambdax]=eig(K(1:45,1:45),M(1:45,1:45));
% [phiy,lambday]=eig(K(46:90,46:90),M(1:45,1:45));
[wnx,d] = spdiags(lambdax);
Fnx=wnx/(2*pi);

in=[100*randn(1,2^16),zeros(1,2^12)];
nt=length(in);
Fs=100;
dt=1/Fs;
t=0:dt:(nt-1)*dt;

fn1in=20*sin(wnx(45)*t);
fn2in=10*sin(wnx(44)*t);
staticin=1000*ones(nt,1);

nenv=2^13;
denv=(pi/2)/nenv;
envec=0:denv:(nenv-1)*denv;
senv=sin(envec);
eenv=cos(envec);
for kk=1:45
invec(kk,:)=in*phix(kk,45)+30*randn(1,nt)+fn1in*phix(kk,45)+fn2in+staticin'*phix(kk,45);
% invec(kk,find(~(in)))=zeros(1,length(find(~(in))));
% end
invec(kk,1:nenv)=senv.*invec(kk,1:nenv);
invec(kk,nt-nenv+1:nt)=eenv.*invec(kk,nt-nenv+1:nt);
end

ninvec=[zeros(45,2^12),invec,zeros(45,2^15)];
nnt=length(ninvec);
tt=0:dt:(nnt-1)*dt;

Xo=zeros(2*n,nnt);
[Ads,Bds]=c2d(As,Bs1,dt);
h=waitbar(0);

for kk=1:nnt-1
    Xo(:,kk+1)=Ads*Xo(:,kk)+Bds*ninvec(:,kk);
    waitbar(kk/(nnt-1),h,'Simulating...');
end
close(h)
ya45=C45*Xo;
da45=cumtrapz(cumtrapz(ya45)*dt)*dt;
%%
close all

gt=0:(nnt-1)*dt;
gvec=1:100:length(tt);
nda45=da45+0.0005*randn(1,length(da45));
figure,plot(tt,ya45,'k'),xlabel('Time (sec)'),ylabel('45F Acceleration (m/sec^2)')
figure,plot(tt,da45,'k'),xlabel('Time (sec)'),ylabel('45F Displacement (m)')
% figure,plot(tt,ninvec',':'),xlabel('Time (sec)'),ylabel('Wind-induced Force (N)')
figure,plot(tt,nda45,':k',gt,da45(gvec),'or'),xlabel('Time (sec)'),ylabel('45F Displacement (m)')
legend('Estimated from Acc.','GPS Measurement')
% figure,plot(tt,

nwdata=load('da45.txt');

figure,plot(tt,da45-nwdata)