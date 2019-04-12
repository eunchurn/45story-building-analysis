clear all; clc; close all;

load patdata;
Fs=3000;
dt=1/Fs;

%% Kalman Filter Design

alpha=0.99999;  % process noise의 비율을 0.00001로 가정한경우
F=[alpha,0,0;1,0,0;0,1,0];  % noise의 Gain이 없는 경우 process noise 무시
H=Fs^2*[1 -2 1];

Q=0.00001; R=0.00001; % process noise, measured noise covariance

y=detrend(patdata{1});

xe=zeros(3,1);
X=zeros(3,length(y));
P(1)=0;

for kk=2:length(y)
  xe=F*X(:,kk-1);
  P=F*P*F'+Q;
  G=P*H'*inv([H*P*H'+R]);
  X(:,kk)=xe+G*(y(kk)-H*xe);
  P=(eye(3)-G*H)*P;
end

%% Ploting States

figure
subplot(211),plot(1:length(X),X'),legend('Priori KF','Kalman Filter Estimated','Posteriori KF')
xlabel('Samples'),ylabel('Estimated Displacement')
subplot(212),plot(y),
xlabel('Samples'),ylabel('Measured Acceleration')