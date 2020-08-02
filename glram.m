% ----- INPUT -----
% X: p * q * n data array
% ptilde: specified dimension of A
% qtilde: specified dimension of B
% ind_mean: 1=centered, 0=non-centered 
%
% ----- OUTPUT -----
% A: extracted left basis (p * ptilde) 
% B: extracted right basis (q * qtilde) 
% s_A: singular values for A
% s_B: singular values for B


function [A,B,s_A,s_B]=glram(X,ptilde,qtilde,ind_mean)

[p,q,n]=size(X);
if ind_mean==1
    mu_X=mean(X,3);
else
    mu_X=0;
end

A1=unifrnd(-1,1,p,ptilde); % initial value
rmsre1=1; d=1; sg=1;

dc=10^(-5);
while d>dc && sg<=30
    A0=A1; rmsre0=rmsre1;
    
    M_B=zeros(q,q);
    for i=1:n
        M_B=M_B+(X(:,:,i)-mu_X)'*(A0*A0')*(X(:,:,i)-mu_X);
    end
    [B1,s_B]=svds(M_B,qtilde);
    
    M_A=zeros(p,p);
    for i=1:n
        M_A=M_A+(X(:,:,i)-mu_X)*(B1*B1')*(X(:,:,i)-mu_X)';
    end
    [A1,s_A]=svds(M_A,ptilde);
    
    for i=1:n
        rmsre1=rmsre1+sum(sum((X(:,:,i)-mu_X...
            -A1*(A1'*(X(:,:,i)-mu_X)*B1)*B1').^2));
    end
    rmsre1=(rmsre1/n)^(0.5);
    d=abs(rmsre0-rmsre1)/rmsre0;
    sg=sg+1;
end
[A,s_A]=svd(M_A/n);  A=A(:,1:ptilde);  s_A=diag(s_A);
[B,s_B]=svd(M_B/n);  B=B(:,1:qtilde);  s_B=diag(s_B);
sg




