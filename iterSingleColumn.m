function [m_hat, r] = iterSingleColumn(A, M, J, k, debug)
    
    [dim, ~] = size(A);
    curr_err = norm(A*M - eye(dim), 'fro');

    E = eye(dim);
    ek = E(:, k);

    % init I and J sets
    % display('I is ');
    I = getIfromJ(A, J);
    [~, n2] = size(J);
    [n1, ~] = size(I);

    % Compute QR of A_hat
    A_hat = A(I, J);
    [Q, R] = qr(A_hat);
    R_hat = R(1:n2, :)
    ek_hat = ek(I);
    
    % update m_hat by the solution of least squares 
    c_hat = Q'*ek_hat;

    m_hat = R_hat\c_hat(1:n2);
    
    % compute residuals
    % display('residual is ');
    r = A(:, J)*m_hat - ek;
    

end
