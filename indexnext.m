function [next] = indexnext(now,n,step)


if step>0
    if now+step>n
        next=indexnext(now-n,n,step);
    else
        next=now+step;
    end
else
    if now+step<1
        next=indexnext(now+n,n,step);
    else
        next=now+step;
    end
end
end

