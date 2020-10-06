function XX = createLinespace(x1,x2,dx)

n = length(dx)+1;
len = sum(dx);
X=zeros(n);
XX = zeros(n-1,3);
X(1) = x1;
for i = 2:n
    X(i)=(x2-x1)/len*dx(i-1)+X(i-1);
    XX(i-1,:)=linspace(X(i-1),X(i),3);
end
XX=XX';

