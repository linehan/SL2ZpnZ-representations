#!/bin/octave-cli

% MODULAR ARITHMETIC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TOTIENT  calculates the totient function (also
%         called the Euler Phi function) of any
%         positive integer n.
%
%Usage:   f = totient(n)
%
%         The totient function computes the number of
%         integers from the set (1...n-1) that are
%         relatively prime to n. It can be used to
%         describe the multiplicative structure of
%         all Galois fields GF(q).
%
%         n can be any size
%
%         Tested under version 5.2.0
%
%Ref: "Error Control Systems" by S.B Wicker
%
%see also Primes, Factor

% Paul Godfrey
% pgodfrey@conexant.com
% 10-23-2001
function [t] = totient(n)
	[r c]=size(n);
	n=reshape(n,1,r*c);
	t=zeros(1,r*c);
	f=zeros(1,10);

	for k=1:r*c;
	        nk=n(k);
	        f=unique(factor(nk));
	        t(k)=nk*prod(1-1./f);
	end

	t=reshape(t,r,c);
	p=find(n==1);
	t(p)=1;
	t=round(t);
	return

	%a demo of this program is
	%n=(1:50).';
	%t=totient(n);
	%[n t]
end

% vpi/powermod: Compute mod(a^d,n)
% usage: R = powermod(a,d,n)
% 
% powermod is MUCH faster than direct exponentiation
% with mod for large numbers. powermod does NOT
% suppoort array or vector inputs, only scalar inputs.
%
% arguments: (input)
%  a,d,n - vpi SCALAR integers, or numeric values
%
% arguments: (output)
%  R - a vpi scalar integer, representing mod(a^d,n)
%
% Example:
%  Compare exponentiation plus mod to
%  the direct application of powermod:
%
%  tic,M = powermod(vpi(123),200,497);toc
%  Elapsed time is 0.044618 seconds.
%
%  tic,M = mod(vpi(123)^200,497);toc
%  Elapsed time is 0.971667 seconds.
%
%
%  See also: power, mod, rem, quotient
%  
% 
%  Author: John D'Errico
%  e-mail: woodchips@rochester.rr.com
%  Release: 1.0
%  Release date: 1/19/09
function R = powermod(a,d,n)
        % convert d to binary, either from a vpi
        % or a double
        if ~isa(d,'vpi')
                db = dec2bin(d);
        else
                db = vpi2bin(d);
        end

        db = fliplr(db == '1');

        % if a is too large, the repeated squarings will
        % cause flint overflow as a double
        if (a > 2^26) || (n > 2^26)
                a = vpi(a);
        end

        if isnumeric(a) && isnumeric(d) && isnumeric(n)
                % pure numeric
                % use the binary expansion of d to form the
                % desired power as efficiently as possible,
                % repeatedly squaring a on each pass.
                if db(1)
                        R = mod(a,n);
                else
                        R = 1;
                end
                for i = 2:length(db)
                        if i > 2
                                a2 = mod(a2*a2,n);
                        else
                                a2 = mod(a*a,n);
                        end
            
                        % do we need to multiply this power
                        % of a into the result?
                        if db(i)
                                % take the mod on each pass through
                                R = mod(R*a2,n);
                        end
                end
        else
                % use the binary expansion of d to form the
                % desired power as efficiently as possible,
                % repeatedly squaring a on each pass.
                if db(1)
                        R = mod(vpi(a),n);
                else
                        R = vpi(1);
                end

                for i = 2:length(db)
                        if i > 2
                                a2 = mod(a2*a2,n);
                        else
                                a2 = mod(a*a,n);
                        end
            
                        % do we need to multiply this power
                        % of a into the result?
                        if db(i)
                                % take the mod on each pass through
                                R = mod(R*a2,n);
                        end
                end
        end
end


function inv = mod_mult_inv(x, q)
	[G, U, V] = gcd(x,q);

	% The inverse of x (mod q) exists only if gcd(x,q) = 1
	if G == 1  
    		inv = mod(U,q);
	else 
		fprintf('No inverse for %d (mod %d)', x, q);
	end
end

function pass = TEST_is_generator(q, g)

	% If g|q, it's definitely not a generator.
	if mod(g, q) == 0
		pass = 0;
		return;
	end

	% The order of the group is equal to the totient function
        for i=2:totient(q)
		if powermod(g,i,q) == g 
			pass = 0;
			return;
		end
        end
	pass = 1;
end


function g = find_multiplicative_generator(q)

	% Loop over candidate generators
	for g=2:q-1
		if TEST_is_generator(q, g) == 1
			return;
		end
	end

	fprintf('no generator for Z/%dZ\n', q);
	g = 0;
end


function x = find_exponent_of_generator(a, g, q)
	for i=1:q-1
		if mod(g^i,q) == a 
			x = i;
			return;
		end
	end

	x = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE OF BASIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
% permutation_from_primitive_root()
% --------------------------------- 
% Generate a permutation matrix corresponding to 
% some given primitive root of a finite field.
%
% @q	: Order of the finite field
% @r	: Primitive q-th root of unity
% Return: (q-1) x (q-1) permutation matrix
% 
% NOTE
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
function M = permutation_from_primitive_root(q, r)

        M = zeros(q-1);

        for i=1:q-1
                M(i, mod(r^i, q)) = 1;
        end
end


%%
% permutation_princ_series_basis()
% -------------------------------- 
% Generate the permutation of the basis elements for
% the vector space Ind(V_{\psi}) as in Rockmore/Lafferty. 
%
% @q	: Order of the finite field
% Return: (q-1) x (q-1) permutation matrix
%
% NOTE
% This changes the order
%
% e_{a^1}, e_{a^2}, e_{a^3}, ..., e_{a^{q-1}}
%
% to
%
% e_{a^2}, e_{a^4}, ..., e_{a^{q-1}}, e_a, e_{a^3}, e_{a^5}, ... e_{a^{q-2}}
%
% Splitting it up into even and odd powers.
%
function P = permutation_princ_series_basis(q)

        P = zeros(q-1);

        n_evens = floor((q-1)/2);

        for i=1:q-1
                if mod(i,2) == 0
                        %i/2
                        P(i/2, i) = 1;
                else
                        ceil(i/2)
                        %n_evens + ceil(i/2)
                        P(n_evens + ceil(i/2), i) = 1;
                end
        end
end



%%
% get_princ_series_cob()
% ----------------------
% Form the change-of-basis matrix for the principal 
% series vector space, relative to some primitive root.
%
% @q	: Order of the field
% @r	: Primitive q-th root of unity
% Return: (q-1) x (q-1) change-of-basis matrix
%
% TODO
% Is the second permutation matrix really independent of the prim root?
function M = get_princ_series_cob(q, r)
        pp = permutation_from_primitive_root(q, r);
        pr = permutation_princ_series_basis(q);

        M = pr * pp;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINCIPAL SERIES REPRESENTATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rep = princ_rep_U(cob, q, a)

        % All representations are linear transformations 
        % in a (q+1)-dimensional vector space. 
        rep = zeros(q + 1);

        % The subgroup U is isomorphic to the finite field Z/qZ.
        for u=1:q
                % Position of u in basis order
                pos_u = u;

                % The shifted value of the 'u'.
                new_u = mod(u-a, q);

                if new_u == 0 
                        % i.e. if 'u' == 'a', then we
                        % send it to the position of
                        % e_0, which is column 'q'.
                        pos_new_u = q;
                else 
                        % Position of the new 'u'.
                        pos_new_u = new_u;
                end

                % We permute pos_u => pos_new_u 
                rep(pos_u, pos_new_u) = 1;
        end

        % e_infty is always stable.
        rep(q+1,q+1) = 1;

        cob(q,q)     = 1;
        cob(q+1,q+1) = 1;

        % Shift into the re-ordered basis relative to a fixed root of unity
        rep = cob * rep * inv(cob);

        return 
end

function rep = princ_rep_T(cob, q, t)

        % All representations are linear transformations 
        % in a (q+1)-dimensional vector space. 
        rep = zeros(q + 1);

        % The subgroup U is isomorphic to the finite field Z/qZ.
        for u=1:q
                pos_u = u;

                % The shifted value of the 'u'.
                new_u = mod((t^2)*u, q);

                if new_u == 0
                        pos_new_u = q;
                else
                        pos_new_u = new_u;
                end

                % We permute pos_u => pos_new_u 
                rep(pos_u, pos_new_u) = 1;
        end

        % e_infty is always stable.
        rep(q+1,q+1) = 1;

        cob(q,q)     = 1;
        cob(q+1,q+1) = 1;

        % Shift into the re-ordered basis relative to a fixed root of unity
        rep = cob * rep * inv(cob);

        return 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETE SERIES REPRESENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n = field_norm(q, t)
	n = powermod(t, q, q);
end

function tr = field_trace(p, n, x)

        tr = x;

        for i=1:(n-1)
                tr = tr + x^i;
        end

        tr = mod(tr, p)
end


function x = chi(j, p, q)
        if p == q
                x = exp(2*pi*1i*j/p);
        else
                x = exp(2*pi*1i*field_trace(j)/p)
        end
end


function j = jay(nu, q, z)

	% Nu is a non-decomposable character

	s = 0;
        qq = q^2;

        c = 0;
	% t is an element of the unit group (F_{q^2})*
	for t=1:qq
                % t needs to be in the unit group, hence can't divide qq.
                %if mod(t, qq) != 0
                        %fprintf('(need %f) field_norm(%d,%d) = %f\n', z, q, t, field_norm(q,t));
                        if field_norm(q, t) == z
                                c = c+1;
                                s = s + chi(powermod(t, q, q)+mod(t, q), q, q) * nu(t, 1);
                        end
                %end
	end

        c

	j = 1/q * s;
end


%
% (WIP) Get the set of nondecomposable characters
% over the unique quadratic field extension of a 
% finite field F_q.
%
function ch = nondecomposable_characters(q)

	% The unique quadratic field extension of F_q is F_{q^2}
	qq = q^2;

	% We will store the characters as column vectors of a matrix.
	ch = zeros(qq);

	% Obtain a generator
	g = find_multiplicative_generator(qq);

	for j=1:qq-1
		for n=1:qq-1
			v = powermod(g, n, qq); % Never zero
                        ch(v, j) = exp(2*pi*1i*j*n*(1/(qq-1)));
		end
	end
        %exit

	% Because F_{q^2} is the unique quadratic extension
	% of F_q, and Gal(F_{q^2}/F_q) is generated by the
	% Frobenius automorphism x->x^q, we know that the unique
	% conjugate of an element in x \in F_{q^2} is given by
	% 
	% 	\overline{x} = x^q.
	% 
	% The conjugate of a character nu of F_{q^2} is defined 
	% pointwise as 
	%
	% 	\overline{nu}(x) = nu(\overline{x})
	%

	%
	% It can be shown that a character nu of F_{q^2} is
	% decomposable if and only if 
	%
	%	nu = \overline{nu}
	%
	% Thus we filter out the decomposable characters:
	
	% ********************
	% FILTER THEM OUT HERE
	% ********************
	c = 0;

	% 0 is always decomposable
	ch(:, qq) = zeros(qq, 1);
	c = c + 1;

	for i=1:qq-1
		% Explain this later. See Reyes paper.
		if mod(q*i, qq-1) == i
			fprintf('%d is decomposable\n', i);
			ch(:, i) = zeros(qq, 1);
			c = c + 1;
		end
	end
	fprintf('removed %d decomposables\n', c);
	
	% Fix this check up.
	(qq-1) - c 
	q^2 - q
end

%ch = nondecomposable_characters(7);

%for i=1:7^2
        %% If the character isn't decomposable
        %if ch(1,i) != 0
                %for x=1:7-1
                        %for y=1:7-1
                                %% mod(x*y,7) is never 0
                                %jay(ch(:,i), 7, mod(x*y, 7))
                        %end
                %end
	%end
%end

%exit


function basis = disc_rep_basis(q, r)
        basis = permutation_from_primitive_root(q,r);
end



function rep = disc_rep_U(cob, q, a)

        rep = zeros(q-1);

        for i=1:q-1
                rep(i, i) = chi(mod(i*a, q), q, q);
        end

        rep = cob * rep * inv(cob);

	% TODO: 
	% The function chi varies? Or not?
end

function rep = disc_rep_T(cob, q, t)

        rep = zeros(q-1);

        % The subgroup U is isomorphic to the finite field Z/qZ.
        for u=1:q-1
                pos_u = u;

                % The shifted value of the 'u'.
                new_u = mod((mod_mult_inv(t,q)^2)*u, q);

                if new_u == 0
                        pos_new_u = q;
                else
                        pos_new_u = new_u;
                end

                % We permute pos_u => pos_new_u 
                rep(pos_u, pos_new_u) = 1;
        end

        % Shift into the re-ordered basis relative to a fixed root of unity
        rep = cob * rep * inv(cob);

	% TODO: 
	% Add the scalar term for each nondecomposable character

        return 
end


function M = disc_rep_W(cob, q, w)

        ch = nondecomposable_characters(q);

        qq = q^2;

        M = zeros(qq);

        for i=1:qq
                % If the character isn't decomposable
                if ch(1,i) != 0
                        for x=1:q-1

                                sub = zeros(q);

                                for y=1:q-1
                                        % Note that mod(x*y,7) is never 0
                                        jj = jay(ch(:,i), q, mod(x*y, q))

                                        %sub(y,y) = 1/ch(y,i) * jj
                                end

                                %M(x:x+q-1, x:x+q-1) = sub;
                        end
                        return;
                end
        end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST SHIT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pass = TEST_reps_U_sum_to_ident(q, r)

        cob = get_princ_series_cob(q, r);

        Z = zeros(q+1);

        % Element of u
        for a=1:q
                Z = Z + princ_rep_U(cob, q, a);
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
        cob = get_princ_series_cob(q,r);

        for t=1:q-1
                RT = princ_rep_T(cob, q, t);
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

function TEST_permute_order(q, r)

        cob = get_princ_series_cob(q, r);

        n_evens = floor((q-1)/2);

        for i=1:q-1
                rr = zeros(q-1,1);

                rr(mod(r^i, q), 1) = 1;

                vect = cob * rr;

                new = find(vect(:,1));

                if mod(i, 2) == 0
                        if new != i/2
                                fprintf('failed\n');
                                pass = 0;
                                return;
                        end
                else
                        if new != n_evens + ceil(i/2)
                                fprintf('failed\n');
                                pass = 0;
                                return;
                        end
                end

                % Should look like order spelled out in paper
                %fprintf('%d => %d => %d\n', i, mod(r^i, q), find(vect(:,1)));
        end

        pass = 1;
        fprintf('passed\n');
        return
end






Q = 7;
R = 3;

A = disc_rep_W(disc_rep_basis(Q,R), Q, 2);
exit


TEST_is_generator(Q, R);
TEST_permute_order(Q, R);
TEST_reps_U_sum_to_ident(Q, R);
TEST_reps_T(Q,R);





A = disc_rep_U(disc_rep_basis(Q,R), Q, 2)

% Lookin' unitary there, good. 
A * ctranspose(A)
ctranspose(A) * A

for i=1:Q-1
	B = disc_rep_T(disc_rep_basis(Q,R), Q, i)
end


exit;

