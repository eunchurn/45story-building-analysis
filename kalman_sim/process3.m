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
time=load(timefile);
accfile=load(accfile);
Fs=100; dt=1/Fs;

A=[0 1;0 0];
B=[0;1];
G=[0;1];

% H=[1,0];
alpha=0.001;
H=[-alpha^2,2*alpha];

Ad=[1 dt;0 1];
Bd=[dt^2/2;dt];
Gd=Bd;

q=0.3;
r=0.3;
Q=[0 0;0 q];
R=r;
Qd=[q*dt^3/3 q*dt^2/2;q*dt^2/2 q*dt];
Rd=r/dt;
xe=zeros(2,1);X=zeros(2,1);P=zeros(2,2);

[K,P,E]=lqe(A,G,H,q,r);
%%
for kk=1:length(meas_a45)-1
    xe(:,kk+1)=Ad*xe(:,kk); %+Bd*meas_a45(kk);
    P=Ad*P*Ad'+Qd;
    K=P*H'*inv(H*P*H'+Rd);
    X(:,kk+1)=xe(:,kk)+K*(meas_a45(kk+1)-H*xe(:,kk));
    P=(eye(2)-K*H)*P;
end
dd=cumtrapz(cumtrapz(meas_a45))/(Fs^2);
figure,plot(1:length(resp_d45),resp_d45,':k',1:length(X),X(1,:),'-r',1:length(xe),xe(1,:),'-k',1:length(dd),dd,'--b')

xe=zeros(2,1);X2=zeros(2,1);P=zeros(2,2);

[K,P,E]=lqe(A,G,H,q,r);
%%
for kk=1:length(meas_a1)-1
    xe(:,kk+1)=Ad*xe(:,kk); %+Bd*meas_a1(kk);
    P=Ad*P*Ad'+Qd;
    K=P*H'*inv(H*P*H'+Rd);
    X2(:,kk+1)=xe(:,kk)+K*(meas_a1(kk+1)-H*xe(:,kk));
    P=(eye(2)-K*H)*P;
end
dd2=cumtrapz(cumtrapz(meas_a1))/(Fs^2);
figure,plot(1:length(resp_d45),resp_d45,':k',1:length(X2),X(1,:)-X2(1,:),'-r',1:length(X2),X2(1,:),'-k',1:length(dd2),dd2,'--b')