% function of signal lasso method using cvx package
% Output is adjacency matrix of network
% Lambda1=0.6;
% Lambda2=0.4;
% Lambda0=0.1;
function [Adj,XSIZE]=sol_cvx_signal_lasso_normal(Fai,y,SIZE,a,Lambda1,Lambda2,Lambda0)
    XSIZE=SIZE;
        cvx_begin quiet
            variable x(SIZE);
            minimize(Lambda0*(Lambda1*norm(x,1)+Lambda2*norm(x-a,1))+square_pos(norm(y-Fai*x,2)));
        cvx_end
    Adj=reshape(x,SIZE,1);
    clear x;
end

