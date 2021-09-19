function sol = ide_solve(idefun,Core,delays_int,history,tspan,options)
    % Function to solve IDEs with Delay Kernel
    solver_name = 'IDE Delay Runge-Kutta';
    sol.solver  = solver_name;

    t0 = tspan(1); tf = tspan(2);
    %======================OPTIONS======================
    Stats      = ideget(options,'Stats','off');
    printstats = strcmp(Stats,'on');
    y0         = ideget(options,'InitialY',history(t0));
    htry       = ideget(options,'InitialStep',[]);
    h          = htry;
    %===================================================
    % Number of equations
    neq = length(y0);
    
    % Stats
    nsteps  = 0;
    %===================================================
    % VIDE Runge-Kutta Tavernini
    A = [0 1/2   0  0
         0   0 1/2  0 
         0   0   0  1     
         0   0   0  0 ];
    
    b = [1/6; 1/3; 1/3; 1/6];
    s = length(b);
    c = [0 1/2 1/2 1];
    d = [1/2; 1/2; 1];
    %===================================================
    nint = length(delays_int(t0)); % Number of integrals
    % Calculate integral (F) in history
    F       = zeros(nint,1);
    
    % Piece from delays(t0) to grid point
    
    % tj, tj_half, tj_1, Core_tj, Core_tj_1, Core_tj_h are always needed to
    % use Simpson's Method for numerical approximations of integrals
    tj      = delays_int(t0);       % Begin
    for ij = 1:nint
        step    = fix((t0 - tj(ij))/h); % The number of memorized intervals of history
        tj_1    = t0 - step*h;      % End
        tj_half = (tj(ij)+tj_1)/2;      % Half segment

        % Calculate kernels at points
        Core_tj      = Core(t0,tj(ij), history(tj(ij)));
        Core_tj_1    = Core(t0,tj_1, history(tj_1));
        Core_tj_h    = Core(t0,tj_half, history(tj_half));

        % Simpson's method
        F(ij)            = F(ij) + int_simpson(tj_1-tj(ij),Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));

        % Main integral over grid points
        Core_tj = Core_tj_1;
        for j = step-1:-1:0
            tj_1         = t0-j*h;
            tj_half      = tj_1 - h/2;

            Core_tj_1    = Core(t0,tj_1, history(tj_1));
            Core_tj_h    = Core(t0,tj_half, history(tj_half));

            F(ij)            = F(ij) + int_simpson(h,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));

            Core_tj      = Core_tj_1;
        end
    end
    %===================================================
    % Initialization | First Step | Y | K
    t(1)   = t0;      sol.x(1)    = t0;
	y(:,1) = y0;      sol.y(:,1)  = y0;
    k      = 1; % Step  

    Y        = zeros(neq,s);
    K        = zeros(neq,s);
    K(:,1,k) = idefun(t(k),y(:,1),F);
    Core_di  = zeros(nint,s);
    %====================================================

	while (t(k) < tf)
        % Last step 
        if (t(k) + h > tf)
            h = tf - t(k);
        end
        %================================================
        Z      = zeros(nint,s);
        Y(:,1) = y(:,k);
        
        % Runge-Kutta steps
        for i = 2 : s
            ti = t(k) + c(i)*h;
            %===================================================
            % Calculate integral (F)
            if i == 2 || i == 4
                F = zeros(nint,1);

                dtk_begin = delays_int(ti); % lower limit of the integral
                for ij = 1:nint
                    if dtk_begin(ij) < t0
                        %===========================================
                        % Integral begins in History
                        step         = fix((t0 - dtk_begin(ij))/htry); % Step of dtk in history

                        % Add piece from dtk_begin to grid point in history
                        tj           = dtk_begin(ij);
                        tj_1         = t0-step*htry;
                        tj_half      = (tj+tj_1)/2;

                        Core_tj      = Core(ti,tj, history(tj));
                        Core_tj_h    = Core(ti,tj_half, history(tj_half));
                        Core_tj_1    = Core(ti,tj_1, history(tj_1));

                        F(ij)            = F(ij) + int_simpson(tj_1-tj,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));

                        % Main integral in history
                        Core_tj = Core_tj_1;
                        for j = step-1:-1:0
                            tj_1      = t0-j*htry;
                            tj_half   = tj_1 - htry/2;

                            Core_tj_h = Core(ti,tj_half, history(tj_half));
                            Core_tj_1 = Core(ti,tj_1, history(tj_1));

                            F(ij)         = F(ij) + int_simpson(htry,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));

                            Core_tj   = Core_tj_1;
                        end

                        % Add integral in solution to t(k)
                        for j = 2:k
                            tj_half   = t(j) - htry/2;
                            y_half    = ntrp3h(tj_half,t(j-1),y(:,j-1),...
                                            K(:,1,j-1),t(j),y(:,j),K(:,1,j));
                            Core_tj_h = Core(ti,tj_half, y_half);
                            Core_tj_1 = Core(ti,t(j), y(:,j));

                            F(ij)         = F(ij) + int_simpson(htry,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));
                            Core_tj   = Core_tj_1;
                        end
                        %===========================================
                    else 
                        %===========================================
                        % Integral only in solution
                        step      = fix((dtk_begin(ij)-t0)/htry + 1); % Step of dtk in solution

                        % Add piece from dtk_begin to grid point in solution
                        tj_half   = (t(step+1) + dtk_begin(ij))/2;
                        y_begin   = ntrp3h(dtk_begin(ij),t(step),y(:,step),...
                                      K(:,1,step),t(step+1),y(:,step+1),K(:,1,step+1));
                        y_begin_h = ntrp3h(tj_half,t(step),y(:,step),...
                                      K(:,1,step),t(step+1),y(:,step+1),K(:,1,step+1));

                        Core_tj   = Core(ti, dtk_begin(ij), y_begin);
                        Core_tj_h = Core(ti,tj_half, y_begin_h);
                        Core_tj_1 = Core(ti, t(step+1), y(:,step+1));

                        F(ij)         = F(ij) + int_simpson(t(step+1)-dtk_begin(ij),Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));

                        % Main integral to t(k)
                        Core_tj = Core_tj_1;
                        for j = step+2:k
                            tj_half   = t(j) - htry/2;
                            y_half    = ntrp3h(tj_half,t(j-1),y(:,j-1),...
                                            K(:,1,j-1),t(j),y(:,j),K(:,1,j));
                            Core_tj_h = Core(ti,tj_half, y_half);
                            Core_tj_1 = Core(ti,t(j), y(:,j));

                            F(ij)         = F(ij) + int_simpson(htry,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));
                            Core_tj   = Core_tj_1;
                        end
                        %===========================================
                    end
                end
                if i == 2
                    F_half = F;
                end
            end
            %===================================================
            if i == 3
                F = F_half;
            end
            
            %Y2-s
            Y(:,i) = y(:,k) + h * ( K(:,1:i-1,k) * A(1:i-1,i) );

            %Z2-s
            
            Core_di(:,i-1)  = Core(t(k)+d(i-1)*h,t(k)+c(i-1)*h,Y(:,i-1));
%             for jj = 1:i-1
%                 Core_di(:,jj)  = Core(t(k)+c(i)*h,t(k)+c(jj)*h,Y(:,jj));
%             end
            Z(:,i)          = h * (Core_di(:,1:i-1) * A(1:i-1,i));

            %K2-s
            K(:, i, k) = idefun(ti, Y(:,i), F+Z(:,i));
        end
        %===============================================
        % Final approximation of RK Method
        t(k+1)       = t(k) + h;
        y(:,k+1)     = y(:,k) + h * (K(:,:,k) * b);
        sol.x(k+1)   = t(k+1);
        sol.y(:,k+1) = y(:,k+1);
        %===============================================   
        % Next step
        % Hermite extrapolation for K(1,k+1)
        y_k_half     = 3/4*y(:,k) + 1/4*y(:,k+1) + h/4*K(:,1,k);
        Core_tj      = Core(t(k+1), t(k),y(:,k));
        Core_tk      = Core(t(k+1), t(k+1),y(:,k+1));
        Core_tk_half = Core(t(k+1), t(k)+h/2, y_k_half);
        F            = F + int_simpson(h,Core_tj,Core_tk,Core_tk_half);
        
        K(:, 1, k+1) = idefun(t(k+1), y(:,k+1), F);
        
        k            = k + 1;
        nsteps       = nsteps + 1;
        %===============================================
    end
    %===================================================
    % Stats
    sol.stats.nsteps = nsteps;
    
    if printstats
        fprintf(getString(message('MATLAB:odefinalize:LogSuccessfulSteps',...
                                                sprintf('%g',nsteps))));
    end
    %===================================================
end

function yint = ntrp3h(tint,t,y,yp,tnew,ynew,ypnew)
    % Hermite extrapolation
    h       = tnew - t;
    s       = (tint - t)/h;
    s2      = s * s;
    s3      = s * s2;
    slope   = (ynew - y)/h;
    c       = 3*slope - 2*yp - ypnew;
    d       = yp + ypnew - 2*slope;
    
    yint    = y + (h*d*s3 + h*c*s2 + h*yp*s);           
end

function f = int_simpson(h, y_begin, y_end, y_half)
    % Simpson's method
    f = h/6*(y_begin + 4*y_half + y_end);
end