using LinearAlgebra

function ode78e(F::Function, t0, tfinal, y0, tol=1e-6, trace=false; evfcn::Function, etol=1e-13)
# [tout, yout, flag] = ode78e(F, t0, tfinal, y0, tol, trace, evfcn, etol)
#
# Runge-Kutta-Fehlberg 7-8 integrator of F with forward and backward
# integration, with optional user specified Event function (evfcn) to
# terminate integration. Directional events supported.
#
# INPUTS:
# F         :function handle for equations of motion
# t0        :initial time  [scalar]
# tfinal    :final time [scalar]
# y0        :initial state [nx1]
# tol       :integrator absolute tolerance [scalar] {1e-6}
# trace     :debugging switch, 1/0 [boolean] {0}
# evfcn     :event function handle {N/A}
# etol      :event function tolerance {1e-13}
#
# OUTPUTS:
# tout      :output time array [mx1]
# yout      :output state array [mxn]
# flag      :true/false indicator of event function tolerance

#
# LOG
# 10/06/1999
# Jonathan Essen
#   Original Code.
# 05/19/2001
# Marc Compere
#   Backward integration.
# 07/20/2015
# Travis Swenson
#   Added event function
#   Improved accuracy and speed of event trigger.
# 08/03/2015
# Brian D. Anderson
#   Added general event function option.
#   Added backward integration capability.
#   Removed superfluous input.
# 06/20/2016
# Brian D. Anderson
#   Added directional event functionality.
#

# initialize continuation variable
integrating = 1;

# The Fehlberg coefficients:
α = [ 2/27 1/9 1/6 5/12 0.5 5/6 1/6 2/3 1/3 1 0 1]'
β = [[2/27       0       0      0        0         0       0         0     0      0     0 0 0 ]
     [1/36       1/12    0      0        0         0       0         0     0      0     0 0 0 ]
     [1/24       0       1/8    0        0         0       0         0     0      0     0 0 0 ]
     [5/12       0       -25/16 25/16    0         0       0         0     0      0     0 0 0 ]
     [.05        0       0      .25      .2        0       0         0     0      0     0 0 0 ]
     [-25/108    0       0      125/108  -65/27    125/54  0         0     0      0     0 0 0 ]
     [31/300     0       0      0        61/225    -2/9    13/900    0     0      0     0 0 0 ]
     [2          0       0      -53/6    704/45    -107/9  67/90     3     0      0     0 0 0 ]
     [-91/108    0       0      23/108   -976/135  311/54  -19/60    17/6  -1/12  0     0 0 0 ]
     [2383/4100  0       0      -341/164 4496/1025 -301/82 2133/4100 45/82 45/164 18/41 0 0 0 ]
     [3/205      0       0      0        0         -6/41   -3/205    -3/41 3/41   6/41  0 0 0 ]
     [-1777/4100 0       0      -341/164 4496/1025 -289/82 2193/4100 51/82 33/164 12/41 0 1 0 ]]'
χ = [ 0 0 0 0 0 34/105 9/35 9/35 9/280 9/280 0 41/840 41/840 ]'
ψ = [ 1 0 0 0 0 0      0    0    0     0     1 -1     -1     ];
pow     = 1 / 8;

# Initializations
t       = t0;
hmax    = (tfinal - t) / 2.5;
# hmin    = (tfinal - t) / 800000000;  # Original Tweak
# hmin    = (tfinal - t) / 8000000000000;  # Travis' Tweak
hmin    = (tfinal - t) / 1e13;  # Brian's Tweak
h       = (tfinal - t) / 100;  # Tweak this?
y       = y0[:];
f       = y * zeros(1,13);
tout    = t;
yout    = y;
τ     = tol * max(norm(y, Inf), 1);
flag    = true;

if trace
    #  clc, t, h, y
    clearconsole(), t, y
end

# evaluate event function at initial state
if event
    [evf,~,dirn]= feval(evfcn,t,y);
end

if tfinal > t0 #FORWARD INTEGRATION
    restep = 0;
    # The main loop
    while (t < tfinal) & (h >= hmin) & integrating
        if t + h > tfinal; h = tfinal - t; end

        # Compute the slopes
        f[:,1]  = F(t,y)
        for j = 1: 12
            f[:, j+1]   = F(t + α[j] * h, y + h * f * β[:, j])
        end

        # Truncation error term
        gamma1  = h * 41 / 840 * f * ψ;

        # Estimate the error and the acceptable error
        δ   = norm(gamma1, Inf);
        τ     = tol*max(norm(y, Inf), 1.0);

        # Update the solution only if the error is acceptable
        if δ <= τ
            t       = t + h;
            y       = y + h * f * χ;
            tout    = [tout; t]
            yout    = [yout; y];
        end
        if trace
            t, y
        end

        # determine event trigger
        restep = 0;
        if (event == 1)
            # compute event function, save previous value
            evf_prev    = evf
            evf, ~, dirn = evfcn(t, y)

            # determine sign change in event function
            cond1   = sign(evf) ~= sign(evf_prev);    #sign switch
            cond2   = sign(evf) == dirn || dirn == 0; #switch direction
            cond3   = sign(evf_prev) ~= 0;            #avoid trigger at prev. exact event
            if cond1 && cond2 && cond3
                # If we need to update algorithm
                if abs(evf) > etol
                    # Delete last step to zoom in and refine
                    tout[end]   = [];
                    yout[end,:] = [];
                    t           = t - h;
                    y           = y - h * f * χ;
                    evf         = evf_prev;
                    restep      = 1;
                else
                    integrating = 0;
                end
            end
        end

        # Update the step size
        if δ ~= 0.0 && restep == 0 #avoid adaptive step size when refining event
            h   = min([hmax, 0.8 * h * (τ / δ) ^ pow]);
        end

        # reduce step size when backtracking to find exact event trigger
        if restep == 1;
            hmax    = h / 5;
            h       = hmax;
        end
    end

    if (event==0)
        if (t < tfinal)
            disp('SINGULARITY LIKELY.')
            t
        end
    else
        if abs(evf) > etol
            # print statement bad for parallel execution
#             disp('EVENT TOLERANCE NOT ACHIEVED.')
#             fprintf('|df|=|#.1e|>#.1e\r\n',evf,etol)
            flag    = false;
        end
    end
else #BACKWARD INTEGRATION
    restep = 0;
    # The main loop
    while (t > tfinal) & (h <= hmin) & integrating
        if t + h < tfinal, h = tfinal - t; end

        # Compute the slopes
        f(:,1)  = feval(F,t,y);
        for j = 1: 12
            f(:, j+1)   = feval(F, t + α(j) * h, y + h * f * β(:, j));
        end

        # Truncation error term
        gamma1  = h * 41 / 840 * f * ψ;

        # Estimate the error and the acceptable error
        δ   = norm(gamma1, 'inf');
        τ     = tol*max(norm(y, 'inf'), 1.0);

        # Update the solution only if the error is acceptable
        if δ <= τ
            t       = t + h;
            y       = y + h * f * χ;
            tout    = [tout; t];
            yout    = [yout; y.'];
        end
        if trace
            t, y
        end

        # determine event trigger
        restep = 0;
        if (event == 1)
            # compute event function, save previous value
            evf_prev    = evf;
            [evf,~,dirn]= feval(evfcn, t, y);

            # determine sign change in event function
            cond1   = sign(evf) ~= sign(evf_prev);    #sign switch
            cond2   = sign(evf) == dirn || dirn == 0; #switch direction
            cond3   = sign(evf_prev) ~= 0;            #avoid trigger at prev. exact event
            if cond1 && cond2 && cond3
                # If we need to update algorithm
                if abs(evf) > etol
                    # Delete last step to zoom in and refine
                    tout(end)   = [];
                    yout(end,:) = [];
                    t           = t - h;
                    y           = y - h * f * χ;
                    evf         = evf_prev;
                    restep      = 1;
                else
                    integrating = 0;
                end
            end
        end

        # Update the step size
        if δ ~= 0.0 && restep == 0 #avoid adaptive step size when refining event
            h   = max([hmax, 0.8 * h * (τ / δ) ^ pow]);
        end

        # reduce step size when backtracking to find exact event trigger
        if restep == 1;
            hmax    = h / 5;
            h       = hmax;
        end
    end

    if (event==0)
        if (t > tfinal)
            disp('SINGULARITY LIKELY.')
            t
        end
    else
        if abs(evf) > etol
            # print statement bad for parallel execution
#             disp('EVENT TOLERANCE NOT ACHIEVED.')
#             fprintf('|df|=|#.1e|>#.1e\r\n',evf,etol)
            flag    = false;
        end
    end
end

return tout, yout, flag

end
