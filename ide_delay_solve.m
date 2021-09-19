function sol = ide_delay_solve(idefun,delays,Core,delays_int,history,tspan,options)
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

IntEqs     = ideget(options,'IntEqs',[]);
%===================================================
% Number of equations
neq = length(y0);

d_t0 = delays(t0,y0);
nz = length(d_t0);

% Stats
nsteps   = 0;
nfevals  = 0;
novsteps = 0;
novinc   = true;
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

%     d_i = delays_int(t0);
%     F = [ -cos(t0)+cos(d_i(1));
%            sin(t0)-sin(d_i(2)) ];
% 	F = cos(t0)-cos(delays_int(t0));
%   F = integral(@(s) cos(s).*exp(t0+s),delays_int(t0),t0);
%d_i = delays_int(t0);
%F = [   integral(@(s) sin(s),d_i(1),t0);
%        integral(@(s) cos(s),d_i(2),t0) ];

% f_i  = @(tx,sx) [ -cos(sx(1));
%                    sin(sx(2)) ];
% f_i  = @(tx,sx) sx*exp(tx);
% f_i  = @(tx,sx) exp(sx + tx)*(cos(sx)/2 + sin(sx)/2);

% f_i  = @(tx,sx) [ exp(2*sx(1))*exp(-tx(1)^2)*exp(sx(1)*tx(1))*exp(tx(1));
%                   exp(-tx(2))*exp(-tx(2)^2)*exp(sx(2)*tx(2))*exp(sx(2)) ];

% f_i  = @(tx,sx) [ -exp(sx(1) - sx(1)*tx(1))/(tx(1) - 1);
%                   -(exp(2*tx(2))*exp(-sx(2)*tx(2))*exp(sx(2)))/(tx(2) - 1) ];
% f_i  = @(tx,sx) -exp(-sx*tx)/tx;
         
% d_i = delays_int(t0);
% F  = f_i(ones(nint,1).*t0,ones(nint,1).*t0) - f_i(ones(nint,1).*t0,d_i);
%===================================================
% Initialization | First Step | Y | K
t(1)   = t0;      sol.x(1)    = t0;
y(:,1) = y0;      sol.y(:,1)  = y0;
k      = 1; % Step  

z = zeros(neq,nz);

for j = 1 : nz
    z(:,j) = history(d_t0(j));
end

Y        = zeros(neq,s);
K        = zeros(neq,s);
K(:,1,k) = idefun(t0,y0,z,F);
nfevals  = nfevals + 1;
Core_di  = zeros(nint,s);
%====================================================

while (t(k) < tf)
    % Last step 
    if (t(k) + h > tf)
        h = tf - t(k);
    end
    %================================================
    Z      = zeros(nint,s);
    Y(:,1,k) = y(:,k);

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

                    F(ij)        = F(ij) + int_simpson(tj_1-tj,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));

                    % Main integral in history
                    Core_tj = Core_tj_1;
                    for j = step-1:-1:0
                        tj_1      = t0-j*htry;
                        tj_half   = tj_1 - htry/2;

                        Core_tj_h = Core(ti,tj_half, history(tj_half));
                        Core_tj_1 = Core(ti,tj_1, history(tj_1));

                        F(ij)     = F(ij) + int_simpson(htry,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));

                        Core_tj   = Core_tj_1;
                    end

%                         F_i = [ -cos(t0)+cos(dtk_begin(1));
%                                   sin(t0)-sin(dtk_begin(2)) ];
%                         F(ij) = F(ij) + F_i(ij);

                    % Add integral in solution to t(k)
                    for j = 2:k
                        tj_half   = t(j) - htry/2;
%                             y_half    = ntrp3h(tj_half,t(j-1),y(:,j-1),...
%                                             K(:,1,j-1),t(j),y(:,j),K(:,1,j))
                        y_begin = zeros(neq,1);
                        y_half  = zeros(neq,1);
                        y_end   = zeros(neq,1);
                        for j1 = 1:neq
                            if ismember(j1,IntEqs)
                                Y_step     = [Y(j1,:,j-1) Y(j1,1,j)];

                                y_begin(j1) = Y_step*W(0);
                                y_half(j1)  = Y_step*W(0.5);
                                y_end(j1)   = Y_step*W(1);
                            else
                                y_begin(j1) = y(j1,j-1);
                                y_half(j1)  = y(j1,j-1) + htry * (K(j1,:,j-1) * b4(1/2));
                                y_end(j1)   = y(j1,j);
                            end
                        end

                        Core_tj   = Core(ti,t(j-1),  y_begin);
                        Core_tj_h = Core(ti,tj_half, y_half);
                        Core_tj_1 = Core(ti,t(j),    y_end);
%                             y_half
                        %y_half = y(:,j-1) + htry * (K(:,:,j-1) * b4(1/2))
                        %y_half
%++++++++++++++++++++++++++++++++++++++++++++++++
%                             y_half    = [ cos(tj_half);
%                                           sin(tj_half) ]
%++++++++++++++++++++++++++++++++++++++++++++++++
                        F(ij)     = F(ij) + int_simpson(htry,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));
%                             for j1 = 1:neq
%                                 if ~ismember(j1,IntEqs)
%                                     Core_tj   = Core_tj_1;
%                                 end
%                             end
                    end

%                         F_i = [ -cos(t(k))+cos(t0);
%                                     sin(t(k))-sin(t0) ];
%                         F(ij) = F(ij) + F_i(ij);
                    %F(ij) = F(ij) -exp(t(k))+exp(t0);
                    %===========================================
                else 
                    %===========================================
                    % Integral only in solution
                    step      = fix((dtk_begin(ij)-t0)/htry + 1); % Step of dtk in solution

                    % Add piece from dtk_begin to grid point in solution
                    tj_half   = (t(step+1) + dtk_begin(ij))/2;

                    y_begin = zeros(neq,1);
                    y_half  = zeros(neq,1);
                    y_end   = zeros(neq,1);
                    for j1 = 1:neq
                        if ismember(j1,IntEqs)
                            Y_step     = [Y(j1,:,step) Y(j1,1,step+1)];

                            y_begin(j1) = Y_step*W((dtk_begin(ij)-t(step))/htry);
                            y_half(j1)  = Y_step*W((tj_half-t(step))/htry);
                            y_end(j1)   = Y_step*W(1);
                        else
                            y_begin(j1) = y(j1,step) + htry * (K(j1,:,step) * b4((dtk_begin(ij)-t(step))/htry));
                            y_half(j1)  = y(j1,step) + htry * (K(j1,:,step) * b4((tj_half-t(step))/htry));
                            y_end(j1)   = y(j1,step+1);
                        end
                    end
%                         y_begin   = y(:,step) + htry * (K(:,:,step) * b4((dtk_begin(ij)-t(step))/htry));
%                         y_begin_h = y(:,step) + htry * (K(:,:,step) * b4((tj_half-t(step))/htry));
%                         y_begin   = ntrp3h(dtk_begin(ij),t(step),y(:,step),...
%                                       K(:,1,step),t(step+1),y(:,step+1),K(:,1,step+1));
%                         y_begin_h = ntrp3h(tj_half,t(step),y(:,step),...
%                                       K(:,1,step),t(step+1),y(:,step+1),K(:,1,step+1));

                    Core_tj   = Core(ti, dtk_begin(ij), y_begin);
                    Core_tj_h = Core(ti, tj_half,       y_half);
                    Core_tj_1 = Core(ti, t(step+1),     y_end);

                    F(ij)     = F(ij) + int_simpson(t(step+1)-dtk_begin(ij),Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));

                    % Main integral to t(k)
                    for j = step+2:k
                        tj_half   = t(j) - htry/2;

                        for j1 = 1:neq
                            if ismember(j1,IntEqs)
                                Y_step     = [Y(j1,:,j-1) Y(j1,1,j)];

                                y_begin(j1) = Y_step*W(0);
                                y_half(j1)  = Y_step*W(0.5);
                                y_end(j1)   = Y_step*W(1);
                            else
                                y_begin(j1) = y(j1,j-1);
                                y_half(j1)  = y(j1,j-1) + htry * (K(j1,:,j-1) * b4(1/2));
                                y_end(j1)   = y(j1,j);
                            end
                        end

                        Core_tj   = Core(ti,t(j-1),  y_begin);
                        Core_tj_h = Core(ti,tj_half, y_half);
                        Core_tj_1 = Core(ti,t(j),    y_end);

                        F(ij)     = F(ij) + int_simpson(htry,Core_tj(ij),Core_tj_1(ij),Core_tj_h(ij));
                    end

%                         d_i = delays_int(ti);
%                         F = [ -cos(t(k))+cos(d_i(1));
%                                sin(t(k))-sin(d_i(2)) ];
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
        Y(:,i,k) = y(:,k) + h * ( K(:,1:i-1,k) * A(1:i-1,i) );
%             for j = 1:neq
%                 if ismember(j,IntEqs)
% %                     f_x     = idefun(ti,Y1(:,j),z,zeros(nint,1));
% %                     ZA = 0;
% %                     for jy = 1:s
% %                         ZA = ZA + A(jy,i) * Z(j,i);
% %                     end
% %                     Y(j,i) = f_x(j) + F(j) + h * ZA;
%                     Y(j,i) = K(j,i-1,k);
%                 else
%                     Y(j,i) = y(j,k) + h * ( K(j,1:i-1,k) * A(1:i-1,i) );
%                 end
%             end
        %Z2-s
%             Core_di(:,i-1)  = Core(t(k)+d(i-1)*h,t(k)+c(i-1)*h,Y(:,i-1,k));
        %INSERT HERE THE RECALCULATION OF CORES FOR THE CURRENT STAGE FOR POUZET
        for jj = 1:i-1
            Core_di(:,jj)  = Core(t(k)+c(i)*h,t(k)+c(jj)*h,Y(:,jj,k));
        end
        Z(:,i)     = h * (Core_di(:,1:i-1) * A(1:i-1,i));

        %Finding delays Z
        d_ti = delays(ti,Y(:,i,k));
        find_z();

%             d_i = delays_int(ti);
%             F = [ -cos(t(k))+cos(d_i(1));
%                    sin(t(k))-sin(d_i(2)) ];
%             if i == 2
%                 F_half = F;
%             end
%         d_i = delays_int(ti);
%         F  = f_i(ones(nint,1).*ti,ones(nint,1).*t(k)) - f_i(ones(nint,1).*ti,d_i);

        %K2-s
        K(:, i, k) = idefun(ti, Y(:,i,k), z, F+Z(:,i));
        nfevals  = nfevals + 1;

        for j = 1:neq
            if ismember(j,IntEqs)
                Y(j,i,k) = K(j,i,k);
            end
        end
    end
    %===============================================
    % Final approximation of RK Method
    t(k+1)       = t(k) + h;
    %y(:,k+1) = y(:,k) + h * (K(:,:,k) * b);
    % Integro-differential equations are calculated first
    y(:,k+1) = zeros(neq,1);
    for j = setdiff(1:neq,IntEqs)
        if ~ismember(j,IntEqs)
            y(j,k+1) = y(j,k) + h * (K(j,:,k) * b);
        end
    end

    % Integral equations are calculated after
    Core_di(:,s)  = Core(t(k+1),t(k)+c(s)*h,Y(:,s,k));
    for j = IntEqs
        if ismember(j,IntEqs)
            f = idefun(t(k+1),y(:,k+1),z,F + h*(Core_di * b));
            nfevals  = nfevals + 1;
            y(j,k+1) = f(j);
%                 y(j,k+1) = Y(j,s);
        end
    end
    %===============================================   
    % Next step
    % Hermite extrapolation for K(1,k+1)
    y_k_half = zeros(neq,1);
    for j = 1:neq
        if ismember(j,IntEqs)
            % Hermite interpolation with Z_left = 0 and right Core_di*b
            % and z_prime Core(...)
            %Z_half = (-cos(t(k)+h/2)+cos(t(k)))/h;

%                 Z_half = 1/4 * h * ((Core_di*b) + Core(t(k),t(k),y(:,k)));
%                 d_ti = delays(t(k)+h/2,Y(:,3,k));
%                 find_z();
%                 f = idefun(t(k)+h/2,y(:,k),z,F_half + Z_half);
%                 y_k_half(j) = f(j);
            y_k_half(j) = [Y(j,:,k) y(j,k+1)] * W(0.5);

            %y_k_half(j) = K(j,3,k);
            %y_k_half(j) = Y(j,3)
            %Y(j,3);
            %cos(t(k)+h/2);
            %y_k_half(j) = cos(t(k)+h/2);
            %y_k_half(j) = exp(t(k)+h/2);
        else
            y_k_half(j)     = y(j,k) + h * (K(j,:,k) * b4(1/2));
            %y_k_half(j) = 3/4*y(j,k) + 1/4*y(j,k+1) + h/4*K(j,1,k);
        end
    end
%         y_k_half;

%         y_k_half = y(:,k) + h * (K(:,:,k) * b4(1/2));
%         y_h(:,k) = y_k_half;

    Core_tj      = Core(t(k+1), t(k),y(:,k));
    Core_tk      = Core(t(k+1), t(k+1),y(:,k+1));
    Core_tk_half = Core(t(k+1), t(k)+h/2, y_k_half);
    %F            = F + int_simpson(h,Core_tj,Core_tk,Core_tk_half);

%         d_i = delays_int(t(k+1));
%         F = [ -cos(t(k+1))+cos(d_i(1));
%                sin(t(k+1))-sin(d_i(2)) ];

%         f1 = int_simpson(h,Core_tj,Core_tk,Core_tk_half);
%         f2 = h*(Core_tj+Core_tk)/2;
%         f3 = [ -cos(t(k+1))+cos(t(k));
%                 sin(t(k+1))-sin(t(k)) ];


    F = F + int_simpson(h,Core_tj,Core_tk,Core_tk_half);
    %F  = F + f_i(ones(nint,1).*t(k+1),ones(nint,1).*t(k+1)) - f_i(ones(nint,1).*t(k+1),ones(nint,1).*t(k));
    
    %F = F + h*(Core_tj+Core_tk)/2;
%         F = F + [ -cos(t(k+1))+cos(t(k));
%                    sin(t(k+1))-sin(t(k)) ];

    d_ti = delays(t(k+1),y(:,k+1));
    find_z();
    K(:,1,k+1) = idefun(t(k+1), y(:,k+1), z, F);
    nfevals  = nfevals + 1;

    for j = 1:neq
        if ismember(j,IntEqs)
            K(j,1,k+1) = y(j,k+1);
        end
    end

    sol.x(k+1)   = t(k+1);
    sol.y(:,k+1) = y(:,k+1);

    k            = k + 1;
    nsteps       = nsteps + 1;
    novinc       = true;
    %===============================================
end
%===================================================
% Stats
sol.stats.nsteps = nsteps;
sol.stats.novsteps = novsteps;

if printstats
    fprintf(getString(message('MATLAB:odefinalize:LogSuccessfulSteps',...
                                            sprintf('%g',nsteps))));
    fprintf(getString(message('MATLAB:odefinalize:LogFunctionEvaluations',... 
                                            sprintf('%g',nfevals))));
end
%===================================================
    function find_z() % Calculate in dt_j
        for kz = 1 : nz
            if d_ti(kz) < t0
                z(:,kz) = history(d_ti(kz));
            elseif ti < d_ti(kz)
                % wrong overlapping
                error('Delays went ahead.');
            elseif t(k) - d_ti(kz) <= 0
                % overlapping
                error('Overlapping.');
            else
                z(:,kz) = find_y_t(d_ti(kz),find_t(d_ti(kz)));
            end
        end
    end

    function y_t = find_y_t(tcur, nstep)
        % find y by given t and step
        theta = (tcur - t(nstep))/htry;
        y_t   = zeros(neq,1);
        for jz = 1:neq
            if ismember(jz,IntEqs)
                Y_step  = [Y(jz,:,nstep) Y(jz,1,nstep+1)];
                y_t(jz) = Y_step*W(theta);
            else
                y_t(jz) = y(jz,nstep) + htry * (K(jz,:,nstep) * b4(theta));
            end
        end
        %y_t = y(:,nstep) + htry * (K(:,:,nstep) * b4(theta));
    end

    function nstep = find_t(tcur)
        % Find step of current t
        %=========Binary search algorithm=========
        iz = 1;
        jz = length(t);
        nstep = fix(jz/2);

        while (((t(nstep+1) < tcur || t(nstep) > tcur)) && iz < jz)
            if tcur > t(nstep)
                iz = nstep + 1;
            else
                jz = nstep - 1;
            end
            nstep = fix((iz+jz)/2);
        end
    end

end

% function yint = ntrp3h(tint,t,y,yp,tnew,ynew,ypnew)
%     % Hermite extrapolation
%     h       = tnew - t;
%     s       = (tint - t)/h;
%     s2      = s * s;
%     s3      = s * s2;
%     slope   = (ynew - y)/h;
%     c       = 3*slope - 2*yp - ypnew;
%     d       = yp + ypnew - 2*slope;
%     
%     yint    = y + (h*d*s3 + h*c*s2 + h*yp*s);           
% end

function f = int_simpson(h, y_begin, y_end, y_half)
    % Simpson's method
    f = h/6*(y_begin + 4*y_half + y_end);
end

function x = b4(a)
    sqrA = a^2;
    x(1,1) = a * (1 + a * (-3/2 + a * 2/3));
    x(2,1) = sqrA * (1 + a * -2/3);
    x(3,1) = sqrA * (1 + a * -2/3);
    x(4,1) = sqrA * (-1/2 + a * 2/3);
end

function x = W(a)
    x(1,1) = 2 * a^2 - 3 * a + 1;
    x(2,1) = a * (-2 * a + 2);
    x(3,1) = x(2,1);
    x(4,1) = -4 * a^2 + 5 * a - 1;
    x(5,1) = 6 * a * (a - 1) + 1;
end