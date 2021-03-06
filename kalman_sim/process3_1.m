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
acc=load(accfile);

Fs=100; dt=1/Fs;
resp_1d=cumtrapz(cumtrapz(acc))/Fs^2;
A=[0 1;0 0];
B=[0;1];
G=[0;1];

% H=[1,0];
alpha=0.001;
H=[-alpha^2,2*alpha];

F=[1 dt;0 1];
Bd=[dt^2/2;dt];
Gd=Bd;

q=0.3;
r=0.3;


Qd=[q*dt^3/3 q*dt^2/2;q*dt^2/2 q*dt];
Rd=r/dt;
x=zeros(2,1);
P=G*q*G';

% [K,P,E1]=lqe(A,G,H,q,r);
%%
for kk=1:length(meas_a45)
    K=P*H'/(H*P*H'+Rd);
    x=x+K*(meas_a45(kk)-H*x);
    P=(eye(2)-K*H)*P;
    ye1(kk)=H*x;
    errcov(kk)=H*P*H';
    
    x=F*x+G*(meas_a45(kk)-meas_a1(kk));
    P=F*P*F'+G*q*G';    
%     xe(:,kk+1)=Ad*xe(:,kk); %+Bd*meas_a45(kk);
%     P=Ad*P*Ad'+Qd;
%     K=P*H'*inv(H*P*H'+Rd);
%     X(:,kk+1)=xe(:,kk)+K*(meas_a45(kk+1)-H*xe(:,kk));
%     P=(eye(2)-K*H)*P;
end
% dd=cumtrapz(cumtrapz(meas_a45))/(Fs^2);
figure,plot(1:length(resp_d45),resp_d45,':k',1:length(ye1),ye1,'-r')

% x=zeros(2,1);

% [K,P,E]=lqe(A,G,H,q,r);
%%
% for kk=1:length(meas_a1)-1
%     K=P*H'/(H*P*H'+Rd);
%     x=x+K*(resp_d1(kk)-H*x);
%     P=(eye(2)-K*H)*P;
%     ye2(kk)=H*x;
%     errcov(kk)=H*P*H';
%     
%     x=F*x+G*meas_a1(kk);
%     P=F*P*F'+G*q*G';   
% end
% % dd2=cumtrapz(cumtrapz(meas_a1))/(Fs^2);
% figure,plot(1:length(resp_d45),resp_d45,':k',1:length(ye2),ye1(1,:)-ye2(1,:),'-r')