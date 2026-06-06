function [y,T]=coupling_dynamic(mysystem,W,y0,T,GAMMA)
% 耦合动力学的模拟过程（龙格库塔）
% 版 权：马闯     %Email：chuangma1990@gmail.com
% %%%%%%%%%%  输入  %%%%%%%%%%%%%%%
% mysystem：节点的动力学系统（函数，mysystem（y））
% W       ：为邻接矩阵中包含耦合强度
% y0      ：系统的初始状态
% T       ：时间序列
% D       ：耦合常数
% GAMMA   ：耦合矩阵
% %%%%%%%%%%  输出  %%%%%%%%%%%%%%%
% y       ：所有节点的在不同时间的状态
%是个三维矩阵，N*M*L（N，节点个数，S，状态的个数，L时间的长度）
%如果存在N，M，L维数为1，三维降为二维
% T       ：时间
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=length(W);
% y0=0.1*ones(n,3)
% delta_t=0.001;
% T=0:delta_t:30;
% D=1
% GAMMA=[0,0,0;
%        0,0,0;
%        0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=RK_4_network(@(t,y0)myfun(mysystem,t,y0,W,GAMMA),T,y0);
y=shiftdim(y);
end

function f=myfun(mysystem,t,y,W,GAMMA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%需要差分的函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W_D = full(W);
    [n,m]=size(y);
    f2=zeros(n,m);
    %%%%%%%%  节点动力学  %%%%%%%%%%%%%%%
    f1=mysystem(t,y);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%   层内耦合   %%%%%%%%%%%%%%%
    for i=1:n
       i_nn=find(W(i,:));%邻居
       for j=i_nn
           f2(i,:)=f2(i,:)-W(i,j)*sum((GAMMA*(y(i,:)-y(j,:))')',1);
       end
    end
    f=f1+f2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function Y=RK_4_network(myfun,T,y0)
%复杂网络上的耦合动力学 龙格库塔 模拟
%版 权：马闯     %Email：chuangma1990@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       输入     %%%%%%%%%%%
%   myfun 自身动力学函数   
%myfun(t,y):t是时间，y是因变量
%   T     时间序列
%   y0    初始每个节点的状态
%%%      输出     %%%%%%%%%%%%
%   y    初始每个节点的状态
%是个三维矩阵，N*M*L（N，节点个数，S，状态的个数，L时间的长度）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [N,M]=size(y0);
    L=length(T);
    Y=zeros(N,M,L);  %初始化
    Y(:,:,1)=y0;
    for i=1:length(T)-1
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        delta_t=T(i+1)-T(i);%时间差
        t=T(i);             %
        y=Y(:,:,i);             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        k1=myfun(t,y);
        k2=myfun(t+delta_t/2,y+delta_t*k1/2);
        k3=myfun(t+delta_t/2,y+delta_t*k2/2);
        k4=myfun(t+delta_t,  y+delta_t*k3);
        Y(:,:,i+1)=Y(:,:,i)+(delta_t/6)*(k1+2*k2+2*k3+k4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

end