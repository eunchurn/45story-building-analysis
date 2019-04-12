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

NFFT=2^nextpow2(length(meas_a45));
dF=Fs/NFFT;
F=0:dF:Fs-dF;
fdata=fft(meas_a45,NFFT);
fdata=fdata-fdata(1);
figure,plot(F,abs(fdata)),xlim([0 F(NFFT/2)])
for kk=2:NFFT/2
%     w2fdata(kk)=fdata(kk)/((2*pi*F(kk))^2);
%     w2fdata(NFFT-kk+2)=fdata(NFFT-kk+2)/((2*pi*F(kk))^2);
    w1fdata(kk)=fdata(kk)/(2*pi*F(kk));
    w1fdata(NFFT-kk+2)=fdata(NFFT-kk+2)/(2*pi*F(kk));
end
% w1fdata=w1fdata-w1fdata(2);
for kk=2:NFFT/2
    w2fdata(kk)=w1fdata(kk)/(2*pi*F(kk));
    w2fdata(NFFT-kk+2)=w1fdata(NFFT-kk+2)/(2*pi*F(kk));
end
% w2fdata=w2fdata(2:end-1);
figure,plot(F,abs(w1fdata)),xlim([0 F(NFFT/2)]),title('Velocity')
figure,plot(F,abs(w2fdata)),xlim([0 F(NFFT/2)]),title('Displacement')
w1data=ifft(w1fdata);
w2data=ifft(w2fdata);

wbias=w2data(length(meas_a45):NFFT);
tb=length(meas_a45):NFFT;
tt=1:length(w2data);
figure,plot(w2data)
figure,plot(tb,wbias)
wff2data=w2data-(fit1(1:length(w2data)))';

figure,plot(1:length(resp_d45),resp_d45,1:length(wff2data),-wff2data)

% 
% for kk=1:NFFT/2
%     xk1(kk)=sin(2*pi*F(kk)*kk*dt);
%     xk2(kk)=cos(2*pi*F(kk+NFFT/2)*kk*dt);
% end
% xk=[xk1 xk2]';
% wk=zeros(NFFT,1);
% mu1=0.1;
% for kk=1:length(meas_a45)
%     e(kk)=fdata(kk)-wk'*xk;
%     wk=wk+2*mu1*xk*e(kk);
% end