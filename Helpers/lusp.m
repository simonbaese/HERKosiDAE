function [ Q,q ] = lusp(A,FV,tol)
%% LU factorization with special pivoting 

% Initialize dimensions.
[n,m] = size(A);

% Initialize indicator for required algebraic components.
fp = flip(find(FV));
r = length(fp);

% Initialize pivot vectors.
p = (1:n)';
q = (1:m)';

% Iteratively change rest matrix.
for k = 1:min(n-1,m)

    % Determine pivot.
    [cv,ri] = max(abs(A(k:n,k:m)));
    
    if k <= r
        ci = fp(k)-k+1;
    else
        [~,ci] = max(cv);
    end
    rp = ri(ci)+k-1;
    cp = ci+k-1;
    
    % Swap rows.
    t = p(k);
    p(k) = p(rp);
    p(rp) = t;
    rt = A(k,:);
    A(k,:) = A(rp,:);
    A(rp,:) = rt;
    
    % Swap columns.
    t = q(k);
    q(k) = q(cp);
    q(cp) = t;
    ct = A(:,k);
    A(:,k) = A(:,cp);
    A(:,cp) = ct;
    
    % Stop factorization if pivot is too small.
    if abs(A(k,k)) >= tol
        rows = (k+1):n;
        cols = (k+1):m;
        A(rows,k) = A(rows,k)/A(k,k);
        A(rows,cols) = A(rows,cols)-A(rows,k)*A(k,cols);
    else
        break
    end
    
end

% Final column swap if m > n.
if n > 0 && m > n
    
    % Determine column pivot.
    if n-1 < r
        ci = fp(n) - n+1;
    else
        [~,ci] = max(abs(A(n,n:m)));
    end
    
    cp = ci + n-1;
    
    % Swap columns.
    t = q(n);
    q(n) = q(cp);
    q(cp) = t;
    
    % Only needed when debugging
    % ct = A(:,n);
    % A(:,n) = A(:,cp);
    % A(:,cp) = ct;
    
end

Q = sparse(q,1:m,1);