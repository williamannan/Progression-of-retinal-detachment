function P_of_s = InterpPoly(s0, X, s)
%% ------- NB: s and s0 must be a row vectors and X must be column(s) vector or matrix  ----------- 

    %%% ------- Inputs-----------------------------
    %%%  s0 - Vector of given Lagrangian or data points
    %%%  X = [X1(s), X2(s)]- matrix of X-coordinates of the data points
    %%%  s - The s-values where the polynomial is to be evaluated

    n = length(s0);          %% Number of points given
    
    % Compute the Lagrange basis polynomials L_i(x) for all i in one go
    L = ones(length(s), n);    %% Matrix to store L_i(x) values
    for i = 1:n
        % Calculate the product for the basis polynomial L_i(x)
        % Avoid division by zero by excluding the current index i
        L(:, i) = prod(bsxfun(@minus, s(:), s0([1:i-1, i+1:end])), 2)./prod(s0(i)-s0([1:i-1, i+1:end]));
    end
    
    % Compute the polynomial values as a weighted sum of basis polynomials
    P_of_s = L*X;
end
