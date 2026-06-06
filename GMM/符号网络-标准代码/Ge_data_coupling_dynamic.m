function [y1,T]=Ge_data_coupling_dynamic(W,S_L,delta_t)
% %%%%%%%%%%%%真实网络
%网络 W
%仿真间隔delta_t
%迭代的步数S_L
%%%%%%%%%%%%%%
n=length(W);
y0=1+2*rand(n,3);
% delta_t=10^-4;
T=0:delta_t:delta_t*S_L;
% D=1;
GAMMA=[1,0,0;
       0,0,0;
       0,0,0];
% [y1,T]=coupling_dynamic(@mysystem2,W,y0,T,D,GAMMA);
[y1,T]=coupling_dynamic(@mysystem2,W,y0,T,GAMMA);

end

function f=mysystem(t,y)
    %节点的动力系统
    %%%%%%%%%%%%%%%%%%%%%%%系统函数
    %Lorenz
    f(:,1)=10*(y(:,2)-y(:,1));
    f(:,2)=28*y(:,1)-(1)*(y(:,1).*y(:,3))-y(:,2);
    f(:,3)=y(:,1).*y(:,2)-(8/3)*y(:,3);
    %%%%%%%%%%%%%%%%%%%%%%
end

function f=mysystem2(t,y)
    %节点的动力系统
    %%%%%%%%%%%%%%%%%%%%%%%系统函数
    %Rossler
    a=0.165;
    b=0.2;
    c=10;
    f(:,1)=-y(:,2)-y(:,3);
    f(:,2)=y(:,1)+a*y(:,2);
    f(:,3)=b+y(:,3).*(y(:,1)-c);
    %%%%%%%%%%%%%%%%%%%%%%
end

