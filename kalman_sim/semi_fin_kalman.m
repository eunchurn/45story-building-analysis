clear all; clc; close all;

measfilea='meas_a45.txt';
measfiled='meas_eq.txt';
respfilea='resp_a45.txt';
respfiled='resp_rd45.txt';
timefile='time.txt';
accfile='elcacc1.txt';

meas_a45=load(measfilea);
meas_a1=load(measfiled);
resp_a45=load(respfilea);
resp_d45=load(respfiled);
t=load(timefile);
accfile=load(accfile);
Fs=100; dt=1/Fs;

alpha=0.001;  % process noise의 비율을 0.00001로 가정한경우
F=[alpha,0,0;1,0,0;0,1,0];  % noise의 Gain이 없는 경우 process noise 무시
G=[1;0;0];
H=Fs^2*[1 -2 1];

q=0.2;r=0.9; % process noise, measured noise covariance

q=0.2;r=0.9;
Qd=[q*dt^5/5 q*dt^4/4 q*dt^3/3;q*dt^4/4 q*dt^3/3 q*dt^2/2;q*dt^3/3 q*dt^2/2 q*dt];
P=G*q*G';
Rd=r/dt;
x=zeros(3,1);

for kk=1:length(meas_a45)

    K=P*H'/(H*P*H'+Rd);
    x=x+K*(resp_d45(kk)-H*x);
    P=(eye(3)-K*H)*P;
    ye(kk)=H*x;
    errcov(kk)=H*P*H';
    
    x=F*x+G*(meas_a45(kk)-meas_a1(kk));
    P=F*P*F'+G*q*G';
end

figure,plot(t,resp_d45,t,ye)
NFFT=2^nextpow2(length(meas_a45));
[Pdxx,F1]=pwelch(resp_d45,[],[],NFFT,Fs);
[Pexx,F2]=pwelch(ye,[],[],NFFT,Fs);

figure,loglog(F1,Pdxx,'-r',F2,Pexx,':k')
figure,plot(errcov)