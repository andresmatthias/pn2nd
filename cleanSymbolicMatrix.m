function C = cleanSymbolicMatrix(CIn)
% CLEANSYMBOLICMATRIX Assuming that entries of a symbolic matrix are 
%   rational functions, find largest coefficient and divide numerator and
%   denominator by this number to avoid expressions of the form 
%   1e35 / (1e35 + 1e35 * sigma_s).
% 
%   In case the entries are not rational functions, one would need to apply
%   some recursive strategy.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

    tolCoeff = 16;
    tolFinal = 16;
    %% extract numerator and denominator
    [n, d] = numden(CIn);
    [~, dn] = numden(n);
    [~, dd] = numden(d);
    C = 0 * CIn;
    sigma_a = str2sym('sigma_a');
    sigma_s = str2sym('sigma_s');
    sigma_t = str2sym('sigma_t');
    if any(dd(:) ~= 1) || any(dn(:) ~= 1)
       error('Matrix entries are not rational functions. Apply recursion from here.') 
    end
    for i = 1 : size(CIn, 1)
        for j = 1 : size(CIn, 2)
            % check, if entry is constant
            if  diff(CIn(i, j), sigma_a) ~= 0 || diff(CIn(i, j), sigma_s) ~=0 || diff(CIn(i, j), sigma_t) ~= 0
                [coefficients_n, terms_n] = coeffs(n(i, j), [sigma_a, sigma_s, sigma_t]); 
                [coefficients_d, terms_d] = coeffs(d(i, j), [sigma_a, sigma_s, sigma_t]); 
                coefficients_n = vpa(coefficients_n, tolCoeff);
                coefficients_d = vpa(coefficients_d, tolCoeff);
                m_n = max(abs(coefficients_n));
                m_d = max(abs(coefficients_d));
                m = max(m_n, m_d);
                coefficients_n = coefficients_n / m;
                coefficients_d = coefficients_d / m;
                C(i, j) = (coefficients_n * terms_n.') / (coefficients_d * terms_d.');
            else
                C(i, j) = CIn(i, j);
            end
        end
    end
    C = vpa(C, tolFinal);
end

