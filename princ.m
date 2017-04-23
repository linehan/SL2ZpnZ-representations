#!/bin/octave-cli

function c = phi(p,j,k)
        c = exp((2*pi*1i*j*k)/(p-1));
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

        n_evens = floor((q-1)/2);

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


function rep = rep_T(cob, q, t)

        % All representations are linear transformations 
        % in a (q+1)-dimensional vector space. 
        rep = zeros(q + 1);

        %for u=1:q
                %fprintf('t=%d, u=%d: %d => %d\n', t, u, mod(u, q), mod((t^2)*u, q));
        %end

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
                        %pos_u = u;%find(cob(u,:));
                        pos_u = find(cob(u,:));
                end

                % The shifted value of the 'u'.
                new_u = mod((t^2)*u, q);

                if new_u == 0 
                        % i.e. if 'u' == 'a', then we
                        % send it to the position of
                        % e_0, which is column 'q'.
                        pos_new_u = q;
                else 
                        % Position of the new 'u'.
                        %pos_new_u = new_u; %find(cob(new_u, :));
                        pos_new_u = find(cob(new_u, :));
                end

                % We permute pos_u => pos_new_u 
                rep(pos_u, pos_new_u) = 1;
        end

        % e_infty is always stable.
        rep(q+1,q+1) = 1;

        %cob(q,q) = 1;
        %cob(q+1,q+1) = 1;

        %rep = cob * rep * inv(cob)

        return 
end

%function rep = rep_T_fast(cob, q, t)
        %% All representations are linear transformations 
        %% in a (q+1)-dimensional vector space. 
        %rep = zeros(q + 1);

        %% The subgroup U is isomorphic to the finite field Z/qZ.
        %for u=1:q
                %pos_u = u;

                %pos_new_u = mod((t^2)*u,q);

                %if pos_new_u == 0
                        %pos_new_u = q;
                %end

                %rep(pos_u, pos_new_u) = 1;
        %end

        %% e_infty is always stable.
        %rep(q+1,q+1) = 1;

        %cob(q,q) = 1;
        %cob(q+1,q+1) = 1;

        %rep = cob * rep;

        %return 
%end


function pass = TEST_reps_U_sum_to_ident(q, r)

        cob = get_cob(q, r);

        Z = zeros(q+1);

        % Element of u
        for a=1:q
                Z = Z + rep_U(cob, q, a);
        end

        if Z(1:q,1:q) == ones(q) && Z(q+1,q+1) == q 
                fprintf('passed\n');
                pass = 1;
        else
                fprintf('failed\n');
                pass = 0;
        end
end


function pass = TEST_reps_T(q, r)
        %
        % PROBLEM: Why are these matrices not of the form given
        % in the paper by Drock? They are passing the test, but
        % still are not in the correct form.
        %
        cob = get_cob(q,r);

        for t=1:q-1
                RT = rep_T(cob, q, t);
                for u=1:q
                        e_u = zeros(q+1,1);
                        e_u(q,1) = 1;
                end

                AUX = RT * e_u;

                need = mod((t^2)*u, q);
                if need == 0
                        need = 7;
                end

                if find(AUX(:,1)) != need
                        fprintf('expected %d, got %d\n', need, find(AUX(:,1)));
                        fprintf('failed\n');
                        pass = 0;
                        return;
                end
        end

        fprintf('passed\n');
        pass = 1;
        return
end
                        
TEST_reps_U_sum_to_ident(7,2);
TEST_reps_T(7,2);

exit;

