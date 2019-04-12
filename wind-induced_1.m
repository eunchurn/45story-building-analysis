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
Bs1=[zeros(n,1);-ones(n/3,1);zeros(n*2/3,1)];
Bs2=[zeros(n,1);-ones(n*2/3,1);zeros(n*1/3,1)];
Cs=[[eye(n,n) zeros(n,n);zeros(n,n) eye(n,n)];As(n+1:2*n,:)];
Ds=zeros(3*n,1);
C45=Cs(2*n+45,:);
D45=0;