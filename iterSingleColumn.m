function [m_hat, pj_sqrs] = iterSingleColumn(A, M, J, k)
    
    [dim, ~] = size(A);
    curr_err = norm(A*M - eye(dim), 'fro');

    E = eye(dim);
    ek = E(:, k);

    % init I and J sets
    I = getIfromJ(A, J);
    [~, n2] = size(J);
    [n1, ~] = size(I);

    % Compute QR of A_hat
    A_hat = A(I, J);
    [Q, R] = qr(A_hat);
    R_hat = R(1:n2, :);
    ek_hat = ek(I);
    
    % update m_hat by the solution of least squares 
    c_hat = Q'*ek_hat;

    m_hat = R_hat\c_hat(1:n2);
    
    % compute residuals
    r = A(:, J)*m_hat - ek;
    norm_r = norm(r)^2;
    
    % compute decrease by taking more things
    pj_sqrs = zeros(dim, 1);
    for j = 1:dim

        pj_sqrs(j) = norm_r - (r'*A(:, j))^2 / norm(A(:, j))^2;

    end

end
