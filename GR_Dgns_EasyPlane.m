clc;
clear all;
close all;

T=10;
epsilon=0.1;
L=10;

mean=zeros(8,1000);
variance=zeros(8,1000);

parfor i=1:8

    for j=1:1:1000

        tic
        j

        if j>1
            X_temp=X(:,10001);
        else
            X_temp=randi([-10,10],2,1);
        end

        X=zeros(2,10001);
        X(:,1)=X_temp;

        for k=1:1:10000

            p=randn(2,1);
            p0=p;
            X(:,k+1)=X(k);
            x=X(:,k);
            for t=1:1:L
                p=p-epsilon/2*[[2*x(1)];[2*x(2)]];
                x=x+epsilon*p;
                p=p-epsilon/2*[[2*x(1)];[2*x(2)]];
            end
            if rand(1,1)<exp(-((x(1)^2+x(2)^2)-(X(1,k)^2+X(2,k)^2))/0.086*T)
                X(:,k+1)=x;
            end


        end

        mean(i,j)=sum(X(1,:))/10001;
        square=0;
        for k=1:1:10001
            square=square+(X(1,k)-mean(i,j))^2;
        end
        variance(i,j)=square/10001;

        toc

    end
end

W=zeros(1,1000);
B=zeros(1,1000);
R=zeros(1,1000);

for i=1:1:1000

    W(i)=sum(variance(:,i))/8;
    square=0;
    mean_B=sum(mean(:,i))/8;
    for j=1:1:8
        square=square+(mean_B-mean(j,i))^2;
    end
    B(i)=square/8*10000;
    R(i)=sqrt((9999/10000*W+1/10000*B)/W);
end


% plot(10000:10000:1000000,R);
% grid on
% ax=gca;
% ax.Title.String='R(EasyPlane,T=10K)';
% ax.Title.FontSize=15;
% ax.Title.FontWeight='Bold';
% ax.XLabel.String='次数';
% ax.YLabel.String='R';
% ax.XLabel.FontSize=12;
% ax.YLabel.FontSize=12;