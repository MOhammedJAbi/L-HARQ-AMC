function [ r ] = bisection_modified( f, a, b, N, eps_step, eps_abs )
    % Check that that neither end-point is a root
    % and if f(a) and f(b) have the same sign, throw an exception.
%     a
        %f(a)
%      b
%f(b)
% 
    if ( abs( f(a) ) < eps_abs )
	r = a;
	return;
    elseif ( abs( f(b) ) < eps_abs )
	r = b;
	return;
    elseif ( f(a)>0 && f(b)>0 )
	r = 0;
	return;
    elseif ( f(a) * f(b) > 0 )
        error( 'f(a) and f(b) do not have opposite signs' );
    end

    % We will iterate N times and if a root was not
    % found after N iterations, an exception will be thrown.

    for k = 1:N
        % Find the mid-point
        c = (a + b)/2;

        % Check if we found a root or whether or not
        % we should continue with:
        %          [a, c] if f(a) and f(c) have opposite signs, or
        %          [c, b] if f(c) and f(b) have opposite signs.

        if ( f(c) == 0 )
            r = c;
            return;
        elseif ( f(c)*f(a) < 0 )
            b = c;
        else
            a = c;
        end

        % If |b - a| < eps_step, check whether or not
        %       |f(a)| < |f(b)| and |f(a)| < eps_abs and return 'a', or
        %       |f(b)| < eps_abs and return 'b'.

        if ( b - a < eps_step )
            if ( abs( f(a) ) < abs( f(b) ) && abs( f(a) ) < eps_abs )
                r = a;
                return;
            elseif ( abs( f(b) ) < eps_abs )
                r = b;
                return;
            end
        end
    end

    error( 'the method did not converge' );
end