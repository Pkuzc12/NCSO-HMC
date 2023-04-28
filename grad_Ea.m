function [gradEa] = grad_Ea(x,n,E_S,field)

gradEa=zeros(4*n*n,1);

S=zeros(6*n*n,1);

grad_S=zeros(6*n*n,1);

grad_field=zeros(6*n*n,1);

trans=zeros(4*n*n,6*n*n);

for i=1:1:n
    for j=1:1:n
        
        index_x=((i-1)*n+j-1)*4;
        index_S=((i-1)*n+j-1)*6;

        theta1=x(index_x+1);
        phi1=x(index_x+2);
        theta2=x(index_x+3);
        phi2=x(index_x+4);

        S(index_S+1:index_S+6)=[sin(theta1)*cos(phi1);sin(theta1)*sin(phi1);cos(theta1);sin(theta2)*cos(phi2);sin(theta2)*sin(phi2);cos(theta2)];
        trans(index_x+1:index_x+2,index_S+1:index_S+3)=[cos(theta1)*cos(phi1),cos(theta1)*sin(phi1),-sin(theta1);-sin(theta1)*sin(phi1),sin(theta1)*cos(phi1),0];
        trans(index_x+3:index_x+4,index_S+4:index_S+6)=[cos(theta2)*cos(phi2),cos(theta2)*sin(phi2),-sin(theta2);-sin(theta2)*sin(phi2),sin(theta2)*cos(phi2),0];
    end
end

for i=1:1:n
    for j=1:1:n

        index=((i-1)*n+j-1)*6;
        grad_field(index+1)=field;
        grad_field(index+4)=field;

    end
end

grad_S=2*E_S*S+grad_field;

gradEa=trans*grad_S;




end

