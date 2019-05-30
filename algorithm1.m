%A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];

A = [ 1.000 0.000 0.000; 5.000 9.000 1000.000; 1.000 2.000 0.000];
[dim, ~] = size(A);

% hardcoded part starts here
% k = 2;
% J = [2 3];
k = 3;
J = [3];
M = eye(dim);
M(3, 2) = 1;
% hardcoded part ends here

[m_hat, pj_sqrs] = iterSingleColumn(A, M, J, k)

% M[J, k] = m_hat
% curr_err = norm(A*M - eye(dim), 'fro')

% running it for the second time! 
%J = [2, J]


pj_sqrs

