%A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];

A = [ 1.000 0.000 0.000; 5.000 9.000 1000.000; 1.000 2.000 0.000];

k = 2;


[dim, ~] = size(A);
M = eye(dim);
M(3, 2) = 1;

curr_err = norm(A*M - eye(dim), 'fro')

E = eye(dim);
ek = E(:, k);

J = [2 3];

% init I and J sets
I = getIfromJ(A, J)
[~, n2] = size(J)
[n1, ~] = size(I)



A_hat = A(I, J);
[Q, R] = qr(A_hat);
R_hat = R(1:n2, :);
ek_hat = ek(I);

c_hat = Q'*ek_hat;


m_hat = R_hat\c_hat(1:n2);

r = A(:, J)*m_hat - ek;
norm_r = norm(r)^2;

pj_sqrs = zeros(dim, 1);
for j = 1:dim
    
    pj_sqrs(j) = norm_r - (r'*A(:, j))^2 / norm(A(:, j))^2;
    
end

pj_sqrs

