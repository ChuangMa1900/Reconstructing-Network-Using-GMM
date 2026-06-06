 function F=clc_F1(w1,w)  
 %w 真实网络
 %w1预测网络
    n=length(w);
    m=sum(sum(w~=0));%真实的边数
    m_n=n*(n-1)-m;%无边的节点对个数
    
    
    tp=length(find(w+w1==2));%预测对的边
    fp=length(find(w-w1==-1));%多预测的边
    
    tpr=tp/m;
    fpr=fp/m_n;
    recall=tpr;  %召回率
    percision=tp/(tp+fp); %精度
    F=2/(1/recall+1/percision);
end