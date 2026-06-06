function F1_err=clc_F1_err_GMM(res_A,QQ,ss,w)
%GMM方法重构网络计算误差与F1
%res_A  %重构邻接矩阵
%QQ     %重构结果
%ss     %重构的耦合强度
%w      %真实网络
%对真实网络进行分层
err=norm(w(:)-res_A(:));

n=length(w);
ss=sort(ss);
w1=zeros(n);
w1(QQ==ss(2))=1;
w1=w1-diag(diag(w1));

w=w-diag(diag(w));
w(w>0)=1;
F=clc_F1(w1,w);
F1_err=[F,err];
