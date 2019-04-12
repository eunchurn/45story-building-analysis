clear all; clc; close all;
% 45층 부산 센텀 리더 건물모델
M=load('m_leader_cen.dat');
C=load('C_leader.dat');
K=load('K_leader.dat');

[r,c]=size(M);
% 자유도
% z = [x1,x2,...,x45,y1,y2,..y45,th1,th2,...,th45]'
% 단위
% ton, meter
n=r;
As=[zeros(n,n) eye(n,n); -inv(M)*K -inv(M)*C];
Bs1=[zeros(n,1);-ones(n/3,1);zeros(n*2/3,1)];
Bs2=[zeros(n,1);-ones(n*2/3,1);zeros(n*1/3,1)];
Cs=[[eye(n,n) zeros(n,n);zeros(n,n) eye(n,n)];As(n+1:2*n,:)];
Ds=zeros(3*n,1);
C45=Cs(2*n+45,:);
C30=Cs(n+30,:);
D45=0;

[phix,lambdax]=eig(K(1:45,1:45),M(1:45,1:45));
[wnx,d] = spdiags(lambdax);
Fnx=wnx/(2*pi);

% NS=load('elcentro_NS.dat');
% EW=load('elcentro_EW.dat');


% NS=[NS(:,1),NS(:,2)*9.8];
% EW=[EW(:,1),EW(:,2)*9.8];

in=load('elcacc1.txt')/1000*9.8;
% in=0.01*randn(10000,1);

nt=length(in);
Fs=100;
dt=1/Fs;
t=0:dt:(nt-1)*dt;

Xo=zeros(2*n,nt);
[Ads,Bds]=c2d(As,Bs1,dt);
for kk=1:nt-1
    Xo(:,kk+1)=Ads*Xo(:,kk)+Bds*in(kk);
end
sys1=ss(As,Bs1,C45,D45);
hz=[0:0.01:50]';
w=hz*2*pi;
% [num,den]=ss2tf(As,Bs,Cs,Ds);
H1=squeeze(freqresp(sys1,w));
figure,semilogy(hz,abs(H1)),grid on
xlabel('Frequency (Hz)'),ylabel('Magnitude')
ya45=C45*Xo;
da45=cumtrapz(cumtrapz(ya45-in')*dt)*dt;

[inr,inc]=size(in);
[ar,ac]=size(ya45);
nin=in+0.001*randn(inr,inc);
nya45=ya45+0.4*randn(ar,ac);
nya30=C30*o+0.4*randn(ar,ac);

nda45=cumtrapz(cumtrapz(nya45-nin')*dt)*dt;
d45=Xo(45,:);

figure,plot(t,in,':k',t,ya45),xlabel('Time(second)'),ylabel('Acceleration (m/sec^2)')
figure,plot(t,da45,':k',t,Xo(45,:)),xlabel('Time(second)'),ylabel('Relative Displacement of 45F (m)')
figure,plot(t,nda45,':k',t,Xo(45,:)),xlabel('Time(second)'),ylabel('Relative Displacement of 45F (m)')
figure,plot(Xo(45,:),Xo(90,:))
ya45=ya45';
d45=d45';

nya45=nya45';
t=t';
% save resp_a45.txt ya45 -ascii
% save resp_rd45.txt d45 -ascii
% save meas_a45.txt nya45 -ascii
% save meas_eq.txt nin -ascii
% save time.txt t -ascii

%% Kalman Filtering

F=Ads;
H=As(n+30,2*n,:);
% [L,P,E] = LQE(A,G,C,Q,R,N);
Q=0.3;R=0.3;N=0;
Rd=R/dt;
Qd=cov(Xo');
P=Bds*Q*Bds';
x=zeros(2*n,1);
%  LQE  Kalman estimator design for continuous-time systems.
%  
%     Given the system
%         .
%         x = Ax + Bu + Gw            {State equation}
%         y = Cx + Du + v             {Measurements}
%  
%     with unbiased process noise w and measurement noise v with 
%     covariances
%  
%         E{ww'} = Q,    E{vv'} = R,    E{wv'} = N ,
%  
%     [L,P,E] = LQE(A,G,C,Q,R,N)  returns the observer gain matrix L
%     such that the stationary Kalman filter
%         .
%         x_e = Ax_e + Bu + L(y - Cx_e - Du)
%  
%     produces an optimal state estimate x_e of x using the sensor
%     measurements y.  The resulting Kalman estimator can be formed
%     with ESTIM.
%  
%     The noise cross-correlation N is set to zero when omitted.  
%     Also returned are the solution P of the associated Riccati 
%     equation
%                              -1
%         AP + PA' - (PC'+G*N)R  (CP+N'*G') + G*Q*G' = 0 
%  
%     and the estimator poles E = EIG(A-L*C).
h=waitbar(0);
for kk=1:length(t)

    K=P*H'/(H*P*H'+Rd);
    x=x+K*([nya45(kk),nya30(kk)-H*x);
    P=(eye(2*n)-K*H)*P;
    ye(kk)=H*x;
    errcov(kk)=H*P*H';
    
    x=F*x+Bds*(nya45(kk)-nin(kk));
    P=F*P*F'+Bds*Q*Bds';
    waitbar(kk/length(t),h)
end
%%
figure,plot(t,d45,t,ye)

nfft=2^nextpow2(length(t));
[Pa,Fa]=pwelch(d45,[],[],nfft,1000);
[Pe,Fe]=pwelch(ye,[],[],nfft,1000);
figure,loglog(Fa,Pa,'-k',Fe,Pe,':r')