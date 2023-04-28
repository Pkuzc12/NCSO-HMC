function [E_S] = grad_Ea_helper(chart,n,parameter)

E_S=zeros(n*n*6);
J1=parameter(1);
J1p=parameter(2);
J2=parameter(3);
J2p=parameter(4);
J3=parameter(5);
K=parameter(6);
G=parameter(7);
Gs=parameter(8);
Kp=parameter(9);
Gp=parameter(10);
Gps=parameter(11);
B=parameter(12);
D=parameter(13);

U=[[1/sqrt(6),-1/sqrt(2),1/sqrt(3)];[1/sqrt(6),1/sqrt(2),1/sqrt(3)];[-sqrt(2/3),0,1/sqrt(3)]];
J=[[0,G,Gs];[G,0,Gs];[Gs,Gs,K]];
Jp=[[0,Gp,Gps];[Gp,0,Gps];[Gps,Gps,Kp]];
order=[[[0,1,0];[0,0,1];[1,0,0]];[[0,0,1];[1,0,0];[0,1,0]];[[1,0,0];[0,1,0];[0,0,1]]];

for i=1:1:n
    for j=1:1:n
        index1=((i-1)*n+j-1)*6;
        E_S(index1+2,index1+2)=B;
        E_S(index1+3,index1+3)=D;
        E_S(index1+5,index1+5)=B;
        E_S(index1+6,index1+6)=D;


        for k=1:1:2
            nei=chart(i,j).left.neighbor(k,:);
            index2=((nei(1)-1)*n+nei(2)-1)*6;
            for t=1:1:3
                E_S(index1+t,index2+t+3)=J1p/2;
            end
            nei=chart(i,j).right.neighbor(k,:);
            index2=((nei(1)-1)*n+nei(2)-1)*6;
            for t=1:1:3
                E_S(index1+t+3,index2+t)=J1p/2;
            end
        end
        nei=chart(i,j).left.neighbor(3,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        for t=1:1:3
            E_S(index1+t,index2+t+3)=J1/2;
        end
        nei=chart(i,j).right.neighbor(3,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        for t=1:1:3
            E_S(index1+t+3,index2+t)=J1/2;
        end


        nei=chart(i,j).left.neighbor2(1,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        for t=1:1:3
            E_S(index1+t,index2+t)=J2p/2;
        end
        
        nei=chart(i,j).right.neighbor2(1,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        for t=1:1:3
            E_S(index1+t+3,index2+t+3)=J2p/2;
        end
        
        nei=chart(i,j).left.neighbor2(4,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        for t=1:1:3
            E_S(index1+t,index2+t)=J2p/2;
        end
        
        nei=chart(i,j).right.neighbor2(4,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        for t=1:1:3
            E_S(index1+t+3,index2+t+3)=J2p/2;
        end


        for k=2:1:3
            nei=chart(i,j).left.neighbor2(k,:);
            index2=((nei(1)-1)*n+nei(2)-1)*6;
            for t=1:1:3
                E_S(index1+t,index2+t)=J2/2;
            end
            nei=chart(i,j).right.neighbor2(k,:);
            index2=((nei(1)-1)*n+nei(2)-1)*6;
            for t=1:1:3
                E_S(index1+t+3,index2+t+3)=J2/2;
            end
        end
        
        for k=5:1:6
            nei=chart(i,j).left.neighbor2(k,:);
            index2=((nei(1)-1)*n+nei(2)-1)*6;
            for t=1:1:3
                E_S(index1+t,index2+t)=J2/2;
            end
            nei=chart(i,j).right.neighbor2(k,:);
            index2=((nei(1)-1)*n+nei(2)-1)*6;
            for t=1:1:3
                E_S(index1+t+3,index2+t+3)=J2/2;
            end
        end



        for k=1:1:3
            nei=chart(i,j).left.neighbor3(k,:);
            index2=((nei(1)-1)*n+nei(2)-1)*6;
            for t=1:1:3
                E_S(index1+t,index2+t+3)=J3/2;
            end
            nei=chart(i,j).right.neighbor3(k,:);
            index2=((nei(1)-1)*n+nei(2)-1)*6;
            for t=1:1:3
                E_S(index1+t+3,index2+t)=J3/2;
            end
        end

       

        nei=chart(i,j).left.neighbor(1,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        E_S(index1+1:index1+3,index2+4:index2+6)=E_S(index1+1:index1+3,index2+4:index2+6)+inv(U)*inv(order(1:3,:))*Jp*order(1:3,:)*U/2;
        nei=chart(i,j).right.neighbor(1,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        E_S(index1+4:index1+6,index2+1:index2+3)=E_S(index1+4:index1+6,index2+1:index2+3)+inv(U)*inv(order(1:3,:))*Jp*order(1:3,:)*U/2;
        nei=chart(i,j).left.neighbor(2,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        E_S(index1+1:index1+3,index2+4:index2+6)=E_S(index1+1:index1+3,index2+4:index2+6)+inv(U)*inv(order(4:6,:))*Jp*order(4:6,:)*U/2;
        nei=chart(i,j).right.neighbor(2,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        E_S(index1+4:index1+6,index2+1:index2+3)=E_S(index1+4:index1+6,index2+1:index2+3)+inv(U)*inv(order(4:6,:))*Jp*order(4:6,:)*U/2;
        nei=chart(i,j).left.neighbor(3,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        E_S(index1+1:index1+3,index2+4:index2+6)=E_S(index1+1:index1+3,index2+4:index2+6)+inv(U)*inv(order(7:9,:))*J*order(7:9,:)*U/2;
        nei=chart(i,j).right.neighbor(3,:);
        index2=((nei(1)-1)*n+nei(2)-1)*6;
        E_S(index1+4:index1+6,index2+1:index2+3)=E_S(index1+4:index1+6,index2+1:index2+3)+inv(U)*inv(order(7:9,:))*J*order(7:9,:)*U/2;


    end
end
end

