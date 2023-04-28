function [Mb] = calculateMb(x,n)

Mb=0;

for i=1:1:n
    for j=1:1:n
        index=((i-1)*n+j-1)*4;
        Mb=Mb+sin(x(index+1))*sin(x(index+2))+sin(x(index+3))*sin(x(index+4));
    end
end

Mb=Mb/(n*n*2);


end

