% function lasso_lasso
clc
clear
load MY
RHO=1;
t_max=1000;
t_max1=1000;

MU=0.001;
MU1=0.01;
MU2=0.1;
LAMBDA=0.1;
T=33


for i=1
    C=M0{i}(1:T,:);
    y=y0{i}(1:T,:);
    [A1,~]=sol_lasso(C,y,RHO,LAMBDA,t_max,MU)

    [B,FitInfo]=lasso(C,y,'cv',5,'Intercept',false,'Lambda',linspace(0,1));
    idxLambda1SE = FitInfo.Index1SE;
    A1 = B(:,idxLambda1SE)



end