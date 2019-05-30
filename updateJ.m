function [J_star] = updateJ(A, r)
    
    [dim, ~] = size(A);
    norm_r = norm(r)^2;
    
    % compute decrease by taking more things
    pj_sqrs = zeros(dim, 1);
    % TODO: Compute J_tilda and iterate over it instead of all the entries
    for j = 1:dim

        pj_sqrs(j) = norm_r - (r'*A(:, j))^2 / norm(A(:, j))^2;

    end
    
    [~, J_star] = min(pj_sqrs)
    