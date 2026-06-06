function F1_err=clc_F1_err_sig_sig(res_A,w,s)
%signal_lasso方法重构网络计算误差与F1
%res_A  %重构邻接矩阵
%w      %真实网络
%s      %真实耦合强度

%对真实网络进行分层
n=length(w);
w_P=zeros(n);
w_N=zeros(n);
w_P(w>0)=1;
w_N(w<0)=1;
%对预测网络进行分层
s=sort(s); %真实值
w1_P=zeros(n);
w1_N=zeros(n);
w1_P(res_A>=(s(2)+s(3))/2)=1;
w1_N(res_A<=(s(2)+s(1))/2)=1;
% %%%计算F1值
F1_P=clc_F1(w1_P,w_P);
F1_N=clc_F1(w1_N,w_N);
F1=(F1_P+F1_N)/2;
F1_err=[F1_P,F1_N,F1,norm(w(:)-res_A(:))];
