#!/bin/octave-cli

function c = phi(p,j,k)
        c = exp((2*pi*1i*j*k)/(p-1));
end


PRIME = 11;
PRIMITIVE_ROOT = 2; % can be anything nonzero above 1 (non-1?)

function xInv = modInv(x,n)
% ModInv(x,n) computes the multiplicative inverse of x modulo n if one
% % exists; errors if no such inverse exists
 if gcd(x,n) ~= 1
     error('x has no inverse modulo n')
end

     [d, a, b]   = gcd(x,n);
     xInv        = mod(a,n);
end




function metadata = make_metadata(prim_root, q)

        metadata = struct(
                'prim_root',    prim_root,
                'q',            q,
                'index_to_value', zeros(1,q), 
                'index_to_power', zeros(1,q), 
                'value_to_index', zeros(1,q), 
                'power_to_index', zeros(1,q),
                'value_to_power', zeros(1,q),
                'mx', zeros(q+1)
        );

        num_distinct_powers = q-1;

        n_evens = floor(num_distinct_powers/2);
        n_odds  = num_distinct_powers - n_evens;

        % pow - power of the primitive root of unity
        % val - value of prim_root^pow (mod q)
        % idx - position of things in the order given in paper 

        for pow=1:num_distinct_powers 

                if mod(pow,2) == 0
                        % Even powers of the root go first, 
                        idx = pow/2;
                else
                        % Then odd powers of the root
                        idx = n_evens + ceil(pow/2);
                end

                val = mod(prim_root^pow, q);

                metadata.index_to_value(idx) = val;
                metadata.index_to_power(idx) = pow;
                metadata.power_to_index(pow) = idx;
                metadata.value_to_index(val) = idx;
                metadata.value_to_power(val) = pow;
        end

        % Choice of primitive root induces a permutation
        % on VALUES 1,...,q-1.
        %R = zeros(q+1);

        %R(q, q)     = 1;
        %R(q+1, q+1) = 1;



        %for pow=1:q+1
                %case i 
                %of 1 do 
                        %% e_0 -> q
                        %metadata.mx(q, i) = 1;
                        %break;
                %of q+1 do
                        %% e_inf -> q+1
                        %metadata.mx(q+1,q+1) = 1;
                        %break;
                %otherwise do
                        %if mod(i,2) = 0
                                %metadata.mx(
                %end
                %if 
                %metadata.mx(i
                 
                
end

%function P = shift_perm_matrix(m, a) 

        %P = zeros(m.q + 1);

        %for i=1:m.q
                %P(mod(i-a, m.q)+1, i) = 1;
        %end
%end


% 
% COLUMNS ARE BASIS "VECTORS" (actually functions)
%
function M = basis_for_induced_vs(meta)

        % The dimensions of the vector space induced in
        % this way is a function of q alone and independent
        % on the choice of character representation \psi.

        basis_matrix = zeros(meta.q + 1);

        % There are q elements of F_q, and one special infty element 
        for pow=1:meta.q-1
                col = meta.power_to_index(pow);
                row = mod(meta.prim_root^pow, meta.q);

                % Write a 1 to the proper position.
                basis_matrix(row, col) = 1;
        end

        % Set a 1 in the position corresponding to e_0.
        basis_matrix(meta.q, meta.q) = 1;

        % Finally, set a 1 in the last position corresponding to e_infty.
        basis_matrix(meta.q+1, meta.q+1) = 1;

        M = basis_matrix;
end


% Remember subgroup U is isomorphic to just the field K as an (additive) group ?
% so u can come from this group

% this should be done as a set of permutation matrices with multiplication
% but until then it is like this...

%
%       Essentially once you fix your initial basis matrix,
%       the shifts by successive u correspond to shifting
%       that matrix UP by one step and wrapping around.
%
%               1 0 0    0 1 0    0 0 0    1 0 0
%               0 1 0 => 0 0 0 => 1 0 0 => 0 1 0
%               0 0 0    1 0 0    0 1 0    0 0 0 
%
function M = representation_U(meta, u)

        rep = zeros(meta.q + 1);

        % There are q elements of F_q, and one special infty element 
        for idx=1:meta.q-1

                v  = meta.index_to_value(idx);
                vv = mod(v - u, meta.q)

                if vv == 0
                        vv = meta.q;
                end

                col = idx;
                row = vv;

                rep(row, col) = 1;
        end
        
        % Finally, set a 1 in the last position corresponding to e_infty.
        rep(meta.q+1, meta.q+1) = 1;

        M = rep;
end


function M = representation_T(meta, a)

        rep = zeros(meta.q + 1);

        % There are q elements of F_q, and one special infty element 
        for idx=1:meta.q-1
        % fix these damn indexes...

                v  = meta.index_to_value(idx)
                vv = mod((a^2)*v, meta.q)

                if vv == 0
                        vv = meta.q; % ??
                end

                ii = meta.value_to_index(vv);

                col = ii;
                row = v;

                %rep(row, col) = phi(modInv(a, meta.q));
                %rep(row, col) = phi(meta.q, 1, -1);
                rep(row,col) = 1;
                %M = phi(modInv(a, meta.q)) * rep;
        end
        
        % Finally, set a 1 in the last position corresponding to e_infty.
        rep(meta.q+1, meta.q+1) = 1; %phi(meta.q, 1, 1);

        M = rep;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% A choice of primitive qth root of unity
% to generate the unit group of the finite
% field is equivalent to a permutation on
% the values 1,...,1-q.
%
% In other words, 
%
%       (Z/qZ)* = <a> 
%
%       a^1, a^2, a^3, ..., a^{q-1}
%       
% We can view 'a' as a map 
%
%       a(k) = a^k
%
% That maps its exponents to values.
%
function P = perm_from_prim_root(r, q)

        P = zeros(q-1);

        for i=1:q-1
                P(mod(r^i, q), i) = 1;
        end
end

function P = permute_reorder(q)

        P = zeros(q-1);

        n_evens = floor((q-1)/2)

        for i=1:q-1
                if mod(i,2) == 0
                        P(i, i/2) = 1;
                else
                        P(i, n_evens + ceil(i/2)) = 1;
                end
        end
end

function M = get_cob(q, r)
        pp = perm_from_prim_root(r, q);
        pr = permute_reorder(q);

        M = pr * pp;
end


function rep = rep_U(cob, q, a)

        % All representations are linear transformations 
        % in a (q+1)-dimensional vector space. 
        rep = zeros(q + 1);

        % The subgroup U is isomorphic to the finite field Z/qZ.
        for u=1:q
                if u == q
                        % The COB matrix is only (q-1)-dimensional,
                        % because it only represents (Z/qZ)*.
                        % 
                        % So we can't ask it for the position of 
                        % u = q. But we don't need to. We interpret 
                        % u = q as u = 0, since q = 0 (mod q), and 
                        % the position of e_0 is fixed at column q 
                        % in our order.
                        pos_u = q;
                else
                        % Position of u in basis order
                        pos_u = find(cob(u,:));
                end

                % The shifted value of the 'u'.
                new_u = mod(u-a, q);

                if new_u == 0 
                        % i.e. if 'u' == 'a', then we
                        % send it to the position of
                        % e_0, which is column 'q'.
                        pos_new_u = q;
                else 
                        % Position of the new 'u'.
                        pos_new_u = find(cob(new_u, :));
                end

                % We permute pos_u => pos_new_u 
                rep(pos_u, pos_new_u) = 1;
        end

        % e_infty is always stable.
        rep(q+1,q+1) = 1;

        return 
end





%for i=1:11-1
        %mod(2^i, 11)
%end


%I = eye(11-1)

%pp = perm_from_prim_root(2, 11)

%pr = permute_reorder(11)

%% this gives us what we have now
%%pp * pr * I

%% So our change of basis matrix is
%cob = pr * pp


%R = zeros(11+1);

%for i=1:11
        %R = R + rep_U(cob, 11, i)
%end


function pass = TEST_reps_U_sum_to_ident(q, r)

        cob = get_cob(q, r);

        Z = zeros(q+1);

        % Element of u
        for a=1:q
                Z = Z + rep_U(cob, q, a);
        end

        if Z(1:q,1:q) == eye(q) & Z(q+1,q+1) == 1 
                pass = 1;
        else
                pass = 0;
        end
end


if TEST_reps_U_sum_to_ident(11, 2) == 1
        fprintf('passed\n');
end



exit

%m = make_metadata(2, 11);


%m.index_to_power
%m.index_to_value



BA = basis_for_induced_vs(m)
P  = shift_perm_matrix(m, 1)

P*BA


%exit

%R = zeros(11+1);

%for i=0:10
        %r = representation_U(m, i)
        %%R = R + r;
%end

%exit

%for i=1:10
        r = representation_T(m, 2)
        %R = R + r;
%end

