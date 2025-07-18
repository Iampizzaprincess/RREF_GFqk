function [R,field] = rref_gfk(A,q,k)
    %% Description:
    % This function converts matrix A to Reduced Row Echelon Form 
    % (RREF) using arithmetic over GF(q^k), where A contains 
    % entries over GF(q^k) and A is in exponential form. This means
    % that if A(1,1) is 1, then this corresponds to the exponent of
    % β^1 where β is a primitive element of GF(q^k). In this 
    % function, -inf corresponds to the 0 element of GF(q^k), since
    % 0 cannot be expressed in exponential form. 
    % 
    % Inputs:
    % A: Matrix in exponential form over GF(q^k)
    % q: base field parameter (should be prime)
    % k: extension field parameter 
    % 
    % Outputs: 
    % R: Matrix in RREF form, expressed exponentially over GF(q^k)
    % field: Outputs a GF(q^k) table. Exponential form can be 
    % mapped to polynomial form using this table. For example,
    % field(i,:) is the polynomial form of β^(i-2) in order of 
    % lowest order to highest order. Note that field(1,:) 
    % corresponds to the 0 element of GF(q^k). 
    %
    % Dependencies:
    % 
    % function gfksub: This function performs a-b over GF(q^k)
    % arithmetic, where a,b are elements GF(q^k). This function is
    % included in this file. 
    %
    % References: 
    % RREF algorithm: Yuanzhen Shao- AN ALGORITHM FOR REDUCING A MATRIX TO
    % ROW ECHELON FORM
    % https://www.math.purdue.edu/~shao92/documents/Algorithm%20REF.pdf
    %
    % Credits: 
    % Author: Jacqueline Malayter
    % Date: July 18, 2025
    %
    % This code was developed by Jacqueline Malayter. If you use or adapt 
    % this code in your own work or publication, please credit the
    % original author by including a reference in your code comments or 
    % documentation.
    %
    % License: This work is licensed under the Creative Commons Attribution
    % 4.0 International License (CC BY 4.0). To view a copy of this 
    % license, visit: https://creativecommons.org/licenses/by/4.0/
    

    % Get primitive polynomial for $\mathbb{F}_{q^{k_i}}$ and
    % create field table
    prim_poly = gfprimfd(k, 'min', q); 
    field = gftuple([-1:q^k-2]',prim_poly,q); % \beta^i corresponds to row i+1 (polynomial form)
    
    % RREF algorithm
    [m, n] = size(A);
    R = A;
    lead = 1;
    for r = 1:m
        if n < lead
            return;
        end
        i = r;
        while R(i, lead) == -inf 
            i = i + 1;
            if i > m
                i = r;
                lead = lead + 1;
                if n < lead
                    return;
                end
            end
        end

        % Swap rows i and r
        temp = R(r, :);
        R(r, :) = R(i, :);
        R(i, :) = temp;

        % Find multiplicative inverse of R(r, lead) mod q^k
        b = 0;
        a = R(r,lead);
        inv = gfdiv(b,a,field);

        % Scale row r to make pivot = 1
        R(r,:) = gfmul(R(r, :), inv * ones(size(R(r,:))), field);

        % Eliminate all other entries in this column
        for i = 1:m
            if i ~= r
                factor = R(i, lead);
                prod1 = gfmul(factor*ones(size(R(r,:))), R(r,:),field);
                for ii = 1 : numel(R(i,:))
                    R(i,ii) = gfksub(R(i,ii), prod1(ii), field, q,k);
                end
            end
        end
        
        lead = lead + 1;
    end
end % end rref_gfqk

function out = gfksub(a,b,field,q,k)
    % Takes two inputs a,b \in GF(q^k) in exponential form, converts
    % to polynomial form using the provided lookup table "field",
    % performs a-b, and returns the output in exponential form.
    % Also requires the base field q and the extension field
    % parameter k. 
    
    % supports scalar inputs only
    assert(isscalar(a));
    assert(isscalar(b));

   % convert a and b to polynomial form
    if b == -inf
        b_poly = zeros([1,k]);
    else
        b_poly = field(b+2,:);
    end

    
    if a == -inf
        a_poly = zeros([1,k]);
    else
        a_poly = field(a+2,:);
    end

    % add a and binv polynomials
    ab = a_poly - b_poly;
    ab = mod(ab,q);

    % convert result to exponential form
    out = find(all(ab==field,2)) - 2;
    if out == -1 
        out = -inf; % zero element matlab representation
    end

end %end gfksub 