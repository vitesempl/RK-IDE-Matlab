function main()
clc;
format long;
%===========================FUNCTIONS======================================
% Output(num,hasDelays,err_calc) - output solution
%             num: number of figure
%        hasDelay: 1 - if equation has discrete delays (z), 0 - if otherwise 
%        err_calc: 1 - do error calculation, 0 - no
%==========================================================================
% head = 'Example 1 (integral+discrete delays)'; disp(head)
% fun = @(t) cos(t);
% 
% tspan = [0 5];
% true_sol = fun(tspan(end));
% 
% idefun = @(t,y,z,int) (1+exp(-pi/2))*y-exp(-pi/2)*z-2*exp(-2*t)*int;
% K = @(t,s,y) y*exp(t+s);
% delays = @(t,y) t-pi/2;
% delays_int = @(t) t-pi/2;
% history = @(t) cos(t);
% 
% options = ideset( 'InitialStep',1e-2,...
%                         'Stats','on'); 
% Output(1,1,0);
%==========================================================================
% head = 'Example 2 (Yukihiko)'; disp(head)
% tspan = [0 15];
% 
% idefun = @(t,y,z,i) [ -2.5*i(1);
%                       -7.5*i(2);
%                        -15*i(3); ];
% K = @(t,s,y) [ sin(y(1));
%                sin(y(2));
%                sin(y(3)) ];
% delays = @(t,y) t-1;
% delays_int = @(t) [ t-1;
%                     t-1;
%                     t-1 ];
% history = @(t) [ 1.5;
%                  1.5;
%                  1.5 ];
% 
% options = ideset( 'InitialStep',1e-2,...
%                         'Stats','off');    
% Output(2,1,0);
%==========================================================================
head = 'Renewal_1'; disp(head)
fun = @(t) [ cos(t);
             sin(t) ];

tspan = [0 5];
true_sol = fun(tspan(end));

idefun = @(t,y,z,i) [ -i(1) + z(2,1);    
                       i(2) + z(2,2) ];
K = @(t,s,y) [ (y(1)^2+y(2)^2)*y(2);
               (y(1)^2+y(2)^2)*y(1) ];
% int((cos(s)^2+sin(s)^2)*cos(s),s)
delays = @(t,y) [ t-pi/2;
                  t-pi   ];
delays_int = @(t) [ t-pi;
                    t-pi/2 ];
history = @(t) [ cos(t)
                 sin(t) ];

options = ideset( 'InitialStep',1e-2,...
                        'Stats','off', ...
                       'IntEqs',1);            
Output(3,1,0);
err_c(4)
%==========================================================================
% head = 'Renewal_1 (reverse)'; disp(head)
% fun = @(t) [ sin(t);
%              cos(t) ];
% 
% tspan = [0 3];
% true_sol = fun(tspan(end));
% 
% idefun = @(t,y,z,i) [  i(2) + z(1,2);
%                       -i(1) + z(1,1) ];
% K = @(t,s,y) [ (y(2)^2+y(1)^2)*y(1);
%                (y(2)^2+y(1)^2)*y(2) ];
% delays = @(t,y) [ t-pi/2;
%                   t-pi   ];
% delays_int = @(t) [ t-pi;
%                     t-pi/2 ];
% history = @(t) [ sin(t)
%                  cos(t) ];
% 
% options = ideset( 'InitialStep',1e-2,...
%                         'Stats','off', ...
%                        'IntEqs',2);       
% Output(4,1,1);
%==========================================================================
head = 'Renewal_2'; disp(head)
fun = @(t) [ exp(t);
             exp(-t) ];

tspan = [0 3];
true_sol = fun(tspan(end));

% i1 = @(t) exp(3*t) - exp(2*t - 2);
% i2 = @(t) 1 - exp(- 2*t - 2);
idefun = @(t,y,z,i) [ i(1) * exp(-t)^2 + exp(t-3) * z(2,1);
                     (i(2) - 1) * exp(t-1)^2 * z(2,3) ];
% idefun = @(t,y,z,i) [ i(1) * exp(-t)^2 + z(1,2) * z(2,1);
%                      (i(2) - 1) * z(1,1)^2 * z(2,3)   ];
% K = @(t,s,y) [ (t+2)*y(1)^2*exp((s-t+1)*t);
%                (t+1)*y(1)^2*y(2)*exp((s-t-1)*t) ];
K = @(t,s,y) [ (t+2)*y(1)^2*exp((s-t+1)*t);
               (t+1)*y(1)^2*exp(-s)*exp((s-t-1)*t) ];
% int((t+2)*exp(s)^2*exp((s-t+1)*t),s)
% int((t+1)*exp(s)^2*exp(-s)*exp((s-t-1)*t),s)
delays = @(t,y) [ t-1;
                  t-3;
                  t-4 ];
delays_int = @(t) [ t-1
                    t-2 ];
history = @(t) [ exp(t)
                 exp(-t) ];

options = ideset( 'InitialStep',1e-2,...
                        'Stats','off', ...
                       'IntEqs',1);            
Output(5,1,0);
err_c(6)
%==========================================================================
% head = 'Renewal_2_ (with dicrete x(tau) )'; disp(head)
% fun = @(t) [ exp(t);
%              exp(-t) ];
% 
% tspan = [0 2];
% true_sol = fun(tspan(end));
% 
% idefun = @(t,y,z,i) [ i(1) * exp(-t)^2 + z(1,2) * z(2,1);
%                      (i(2) - 1) * z(1,1)^2 * z(2,3) ];
% K = @(t,s,y) [ (t+2)*y(1)^2*exp((s-t+1)*t);
%                (t+1)*y(1)^2*y(2)*exp((s-t-1)*t) ];
% delays = @(t,y) [ t-1;
%                   t-3;
%                   t-4 ];
% delays_int = @(t) [ t-1
%                     t-2 ];
% history = @(t) [ exp(t)
%                  exp(-t) ];
% 
% options = ideset( 'InitialStep',1e-2,...
%                         'Stats','off', ...
%                        'IntEqs',1);            
% Output(6,1,1);
%==========================================================================
% head = 'Renewal_2 (with SD)'; disp(head)
% fun = @(t) [ exp(t);
%              exp(-t) ];
% 
% tspan = [0 5];
% true_sol = fun(tspan(end));
% 
% idefun = @(t,y,z,i) [ i(1) * exp(-t)^2 + exp(t-3) * z(2,1);
%                      (i(2) - 1) * exp(t-1)^2 * z(2,3) ];
% K = @(t,s,y) [ (t+2)*y(1)^2*exp((s-t+1)*t);
%                (t+1)*y(1)^2*y(2)*exp((s-t-1)*t) ];
% delays = @(t,y) [ t-1;
%                   t-3;
%                   t-4 ];
% delays_int = @(t) [ t-1
%                     t-2 ];
% history = @(t) [ exp(t)
%                  exp(-t) ];
% 
% options = ideset( 'InitialStep',5*1e-3,...
%                         'Stats','off', ...
%                        'IntEqs',1);            
% Output(7,1,1);
%==========================================================================
% head = 'Renewal_0'; disp(head)
% 
% tspan = [0 5];
% 
% idefun = @(t,y,z,i) [ -i(1) + z(2,1) - i(2);    
%                        i(3) + z(2,2) + 2*i(4) ];
% K = @(t,s,y) [ cos(s)*sin(s)*y(1);
%                (1+2*sin(s)^2)*y(2); 
%                (1+2*cos(s)^2)*y(1); 
%                cos(s)*sin(s)*y(2) ];
% % int((cos(s)^2+sin(s)^2)*cos(s),s)
% delays = @(t,y) [ t-pi/2;
%                   t-pi   ];
% delays_int = @(t) [ t-pi;
%                     t-pi;
%                     t-pi/2;
%                     t-pi/2];
% history = @(t) [ cos(t)
%                  sin(t) ];
% 
% options = ideset( 'InitialStep',1e-3,...
%                         'Stats','off', ...
%                        'IntEqs',1);            
% Output(8,1,0);
%err_c(3)
%==========================================================================
% head = 'IE'; disp(head)
% fun = @(t) exp(t);
% 
% tspan = [0 5];
% true_sol = fun(tspan(end));
% 
% idefun = @(t,y,z,i) i;
% K = @(t,s,y) y^2 * exp(t-2*s);
% % int(exp(s)^2 * exp(t-2*s),s)
% delays = @(t,y) t-2;
% delays_int = @(t) t-1;
% history = @(t) exp(t);
% 
% options = ideset( 'InitialStep',1e-2,...
%                         'Stats','off', ...
%                        'IntEqs',1);            
% Output(8,1,1);
%==========================================================================
% head = '2 integrals'; disp(head)
% fun = @(t) exp(t);
% 
% tspan = [0 5];
% true_sol = fun(tspan(end));
% 
% idefun = @(t,y,z,int) exp(1)-exp(t^2)/(z^2)*(int(1)-exp(-2*t)*int(2))*(t-1);
% K = @(t,s,y) [  y*exp(-s*t);
%                 y*exp(t*(2-s))];
% % f_i  = @(tx,sx) [ -exp(sx(1) - sx(1)*tx(1))/(tx(1) - 1);
% %                   -(exp(2*tx(2))*exp(-sx(2)*tx(2))*exp(sx(2)))/(tx(2) - 1) ];
% delays = @(t,y) t-1;
% delays_int = @(t) [ t-1;
%                     t-2];
% history = @(t) exp(t);
% 
% % options = ideset( 'InitialStep',1e-3,...
% %                         'Stats','on');            
% Output(9,1,0);
%==========================================================================
% head = '2 discrete delays + State-Dependent'; disp(head)
% fun = @(t) exp(-t);
% 
% tspan = [0 5];
% true_sol = fun(tspan(end));
% 
% idefun = @(t,y,z,int) -z(1)^((t+1)/2) * exp(-(t-1) / 4) * y^2 * exp(1/2*t^2 + t - 1/2);
% K = @(t,s,y) y*exp(s-s*t);
% % f_i  = @(tx,sx) -exp(-sx*tx)/tx;
% delays = @(t,y) [ (log(y))^2 / (t+1) - 1/2;
%                   (t-1) / 4 ];
% delays_int = @(t) t/2-1;
% history = @(t) exp(-t);
% 
% % options = ideset( 'InitialStep',1e-3,...
% %                         'Stats','on');            
% Output(10,1,0);
% % err_c(11);
%==========================================================================

function Output(num,hasDelays,err_calc)
    if hasDelays
        sol = ide_delay_solve(idefun,delays,K,delays_int,history,tspan,options);
    else
        sol = ide_solve(idefun,K,delays_int,history,tspan,options);
    end
    nz = length(history(tspan(1)));
    if nz == 1
        fprintf('y(%d) = %.8f\n',tspan(2),sol.y(:,end))
    else
        for i = 1:nz
            LegendsStrings{i} = ['y',num2str(i)];
            fprintf('y%d(%d) = %.8f\n',i,tspan(2),sol.y(i,end))
        end
    end
    fprintf('\n')
%     fprintf('Error:')
%     abs(true_sol - sol.y(end))

    close(figure(num))
    f = figure(num);
    pos = get(f,'Position');
    % [left bottom width height]
    if err_calc
        if ispc
            pos = pos + [0 -pos(2)*0.75 pos(3)*1.5 pos(4)*0.5];
        else
            pos = pos + [0 0 0 pos(4)*1.5];
        end
        subplot(1,2,1);
    else
        if ispc
            pos = pos + [0 -pos(2)*0.75 0 pos(4)*0.2];
        else 
            pos = pos + [0 0 0 pos(4)*0.2];
        end
    end
    set(f,'Position',pos);

    plot(sol.x,sol.y,'LineWidth',3);
    %title(head, 'FontSize', 14);
    grid on;
    xlabel('TIME','FontSize',12); ylabel('SOLUTION','FontSize',12);

    if nz ~= 1
        legend(LegendsStrings)
        legend({'x(t)';'y(t)'},'FontSize',16)
        %legend({'r=2.5';'r=7.5';'r=15'},'FontSize',16,'Location','southwest')
    end

    if err_calc
        err = [];
        nsteps = [];

        nstart = 3
        n = 9;
        for steppow = nstart:n
            step = 2^(-steppow);
            options = ideset(options,'InitialStep',step);

            if hasDelays == 1
                sol = ide_delay_solve(idefun,delays,K,delays_int,history,tspan,options);
            else
                sol = ide_solve(idefun,K,delays_int,history,tspan,options);
            end

            err    = [err,abs(true_sol - sol.y(:,end))];
            nsteps = [nsteps,step];

            %abs(true_sol - sol.y(end))
            if steppow > nstart
                %(log10(err(end))-log10(err(end-1)))/(log10(2^(-steppow))-log10(2^(-steppow+1)))
            end
        end

        fprintf('Convergence order: %.4f\n\n',(log10(err(:,end))-log10(err(:,end-1)))/(log10(2^(-n))-log10(2^(-n+1))));

        subplot(1,2,2);
%            loglog(rtol,ourerr,'k','LineWidth',2); 
        loglog(nsteps,err,'k','LineWidth',2); 
        xlabel('STEPSIZE'); ylabel('ERROR');

        ax = gca;
        ax.YAxis.MinorTickValues = kron(logspace(-2,4,7),(4:3:10));

        grid on;
        set(gca,'FontSize',12);
    end
end

function err_c(num)
    err = [];
    nsteps = [];
    
    nstart = 4;
    n = 9;
    for steppow = nstart:n
        step = 2^(-steppow);
        options = ideset(options,'InitialStep',step);

        sol = ide_delay_solve(idefun,delays,K,delays_int,history,tspan,options);

        err    = [err,abs(true_sol - sol.y(:,end))];
        nsteps = [nsteps,step];

        (tspan(2)-tspan(1))/2^-(steppow)
        abs(true_sol - sol.y(:,end))
        if steppow > nstart
            (log10(err(:,end))-log10(err(:,end-1)))/(log10(2^(-steppow))-log10(2^(-steppow+1)))
        end
    end

    fprintf('Convergence order: %.4f\n\n',(log10(err(:,end))-log10(err(:,end-1)))/(log10(2^(-n))-log10(2^(-n+1))));
    close(figure(num));
    f = figure(num);
    
    pos = get(f,'Position');
    % [left bottom width height]
    if ispc
        pos = pos + [0 -pos(2)*0.75 0 pos(4)*0.2];
    else 
        pos = pos + [0 0 0 pos(4)*0.2];
    end
    set(f,'Position',pos);
    %plot(log10(nsteps),log10(err),'k','LineWidth',3)
    loglog(nsteps,err,'k','LineWidth',3); 

    ax = gca;
    ax.YAxis.MinorTickValues = kron(logspace(-2,4,7),(4:3:10));

    grid on;
    xlabel('STEPSIZE'); ylabel('ERROR');
    set(gca,'FontSize',12);
end

end

