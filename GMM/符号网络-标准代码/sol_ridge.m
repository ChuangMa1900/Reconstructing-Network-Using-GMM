function [A_i,res]=sol_ridge(C,y,LAMBDA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%求解带有稀疏约束的线性方程组：min ||C*A_i-y||^2_2+LAMBDA|A_i|
%Copyright：马闯 chuang_m@126.com
%输入：C,y，其中CA_i=y
%      RHO,惩罚参数，一般取1
%      LAMBDA，一范数约束系数
%      t_max，最大迭代次数    
%      MU，   收敛误差      
%输出：A_i：方程组的解
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RHO=1;
% MU=0.001;
% LAMBDA=0.01;
% t_max=100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(C);%n表示未知量个数

A_i=(C'*C+LAMBDA*eye(n))\(C'*y);
res=0;
% A_i=Z;%把Z作为输出结果
end