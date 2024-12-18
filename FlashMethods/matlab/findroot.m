function[X] = findroot(z,K)
 

%s = 1./(1-K);
Fv_min=1/(1-max(K));
Fv_max=1/(1-min(K));
%a=-0.00001;
%b=1.0000001;

a=Fv_min+0.00001;
b=Fv_max-0.00001;


n=size(z);
X=(a+b)/2;

fa=sum(z.*(K-1)./(1+a*(K-1)));
fb=sum(z.*(K-1)./(1+b*(K-1)));
fX=sum(z.*(K-1)./(1+X*(K-1)));
% sum(fa)
% sum(fb)
% sum(fX)

while (abs(a-b)>0.000001)
fa=sum(z.*(K-1)./(1+a*(K-1)));
fb=sum(z.*(K-1)./(1+b*(K-1)));
fX=sum(z.*(K-1)./(1+X*(K-1)));
if (fa*fX<0)
    b=X;
else
    if(fb*fX<0)
        a=X;
    end
end
X=(a+b)/2;
end

