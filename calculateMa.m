function [Ma] = calculateMa(x,n)

Ma=0;

for i=1:1:n
    for j=1:1:n
        index=((i-1)*n+j-1)*4;
        Ma=Ma+sin(x(index+1))*cos(x(index+2))+sin(x(index+3))*cos(x(index+4));
    end
end

Ma=Ma/(n*n*2);


end

