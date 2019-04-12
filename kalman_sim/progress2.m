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

alpha=0.9;  % process noise의 비율을 0.00001로 가정한경우
F=[alpha,0,0;1,0,0;0,1,0];  % noise의 Gain이 없는 경우 process noise 무시
G=[1;0;0];
H=Fs^2*[1 -2 1];

Q=0.02; R=0.02; % process noise, measured noise covariance
% R=eye(3);
S=zeros(3);
% R=[0.02 0 0;0 0.02 0;0 0 0.02];
% S=zeros(3,3);
S=0.2;
% ye=zeros(3,1);
X=zeros(3,length(meas_a45));
% P(1)=0;
P=zeros(3,3);
q=0.3;r=0.5;
% Q=[0 0;0 q];
% Qd=[q*dt^5/5 q*dt^4/4 q*dt^3/3;q*dt^4/4 q*dt^3/3 q*dt^2/2;q*dt^3/3 q*dt^2/2 q*dt];
% P=G*q*G';
Rd=r/dt;
xe=zeros(3,1);

for kk=1:length(meas_a45)
%     xe(:,kk+1)=F*xe(:,kk);
    e=meas_a45(kk)-H*xe(:,kk);
    R=H*P*H'+R;
%     K=(F*P*H'+G*S)*R; 
    K=P*H'/(H*P*H'+Rd);
    xe(:,kk+1)=F*xe(:,kk)+K*e;
    P=F*P*F'+G*Q*G'-K*R*K';
%     x=x+K*(resp_d45(kk)-H*x);
%     P=(eye(3)-K*H)*P;
%     ye(kk)=H*x;
%     errcov(kk)=H*P*H';
    
%     x=F*x+G*meas_a45(kk);
%     P=F*P*F'+G*q*G';
end
%     R=H*P*H'+R;
%     K=(F*P*H'+G*S)*R;
%     xe(:,kk+1)=F*xe(:,kk)+K*e;
%     P=F*P*F'+G*Q*G'-K*R*K';
%     xe=F*X(:,kk-1);
%     P=F*P*F'+Q;
%     G=P*H'*inv([H*P*H'+R]);
%     X(:,kk)=xe+G*(meas_a45(kk)-H*xe);
%     P=(eye(3)-G*H)*P;
% end
figure,plot(t,resp_d45,t,xe(1,2:end))
% 
% clear P;
% clear G;
% xeb=zeros(3,1);
% Xb=zeros(3,length(meas_a1));
% P(1)=0;
% 
% for kk=2:length(meas_a1)
%     xeb=F*Xb(:,kk-1);
%     P=F*P*F'+Q;
%     G=P*H'*inv([H*P*H'+R]);
%     Xb(:,kk)=xeb+G*(meas_a1(kk)-H*xeb);
%     P=(eye(3)-G*H)*P;
% end
% 
% figure,plot(1:length(X),(X-Xb)',1:length(resp_d45),resp_d45),legend('Kalman Filtered Gain1','Kalman Filtered Gain2','Kalman Filtered Gain3','Numerical Data')
% figure,plot(1:length(meas_a45),cumtrapz(cumtrapz(meas_a45)*dt)*dt)
