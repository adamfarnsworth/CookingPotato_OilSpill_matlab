clear;

% domain
xL = 0;
xR = 1;

% diffusion coefficient
D = 0.6;

% boundary conditions
cL = @(t) 0;
cR = @(t) 0;

% source term
f = @(t,x) 0;

% initial conditions
x0 = 0.5;
tau = 0.001;
c_start = @(x) sqrt(4*pi*D*tau)/sqrt(4*pi*D*(tau)) * exp(-(x-x0)^2/(4*D*(tau)));

% time interval
t_start = 0;
t_final = 0.005;

% exact solution
tau = 0.001;
c_exact = @(t,x) sqrt(4*pi*D*tau)/sqrt(4*pi*D*(t+tau)) * exp(-(x-x0)^2/(4*D*(t+tau)));

% space discretization
Nx = 25;

num_splits = 6;

for k = 1:num_splits
    
    x = linspace(xL, xR, Nx);
    dx = (xR-xL)/(Nx-1);
    
    % time-step from stability
    dt = 0.01*dx;
    
    c_im_old = zeros(1,Nx);
    c_im_new = zeros(1,Nx);
    c_cn_old = zeros(1,Nx);
    c_cn_new = zeros(1,Nx);
    c_exa = zeros(1,Nx);
    
    for i = 1:Nx
        c_im_old(i) = c_start(x(i));
        c_cn_old(i) = c_start(x(i));
    end
    
    t = t_start;
    
    A_im = sparse(Nx,Nx);
    A_cn = sparse(Nx,Nx);
    RHS_im = zeros(Nx,1);
    RHS_cn = zeros(Nx,1);
    
    % Crank-Nicolson scheme
    % internal points
    for i = 2:Nx-1
        A_cn(i,i) = 1+dt*D/dx/dx;
        A_cn(i,i-1) = -.5*dt*D/dx/dx;
        A_cn(i,i+1) = -.5*dt*D/dx/dx;
    end
    
    % boundary points
    A_cn(1,1) = 1;
    A_cn(Nx,Nx) = 1;
    
    % fully-implicit scheme
    % internal points
    for i = 2:Nx-1
        A_im(i,i) = 1+2*dt*D/dx/dx;
        A_im(i,i-1) = -dt*D/dx/dx;
        A_im(i,i+1) = -dt*D/dx/dx;
    end
    
    % boundary points
    A_im(1,1) = 1;
    A_im(Nx,Nx) = 1;
    
    while t < t_final
        
        if t + dt > t_final
            dt = t_final-t;
            
            % Crank-Nicolson scheme
            % internal points
            for i = 2:Nx-1
                A_cn(i,i) = 1+dt*D/dx/dx;
                A_cn(i,i-1) = -.5*dt*D/dx/dx;
                A_cn(i,i+1) = -.5*dt*D/dx/dx;
            end
            
            % fully-implicit scheme
            % internal points
            for i = 2:Nx-1
                A_im(i,i) = 1+2*dt*D/dx/dx;
                A_im(i,i-1) = -dt*D/dx/dx;
                A_im(i,i+1) = -dt*D/dx/dx;
            end
        end
        
        % Crank-Nicolson scheme
        % internal points
        for i = 2:Nx-1            
            RHS_cn(i) = c_cn_old(i) + 0.5*dt*(f(t,x(i))+f(t+dt,x(i))) + .5*dt*D*(c_cn_old(i-1)-2*c_cn_old(i)+c_cn_old(i+1))/dx/dx;
        end
        
        % boundary points
        RHS_cn(1) = cL(t+dt);
        RHS_cn(Nx) = cR(t+dt);
        
        % solve system of equations
        c_cn_new = A_cn\RHS_cn;
        c_cn_old = c_cn_new;
        
        % fully-implicit scheme
        % internal points
        for i = 2:Nx-1            
            RHS_im(i) = c_im_old(i) + dt*f(t+dt,x(i));
        end
        
        % boundary points
        RHS_im(1) = cL(t+dt);
        RHS_im(Nx) = cR(t+dt);
        
        % solve system of equations
        c_im_new = A_im\RHS_im;
        c_im_old = c_im_new;
        
        t = t+dt;
        
%         plot(x,c_exa,'LineWidth',2);
%         hold on
%         plot(x,c_cn_new,'o-','LineWidth',1);
%         plot(x,c_im_new,'*-','LineWidth',1);
%         hold off
%         xlabel('x');
%         ylabel('c');
%         axis([xL xR 0 1]);
%         pause(1000*dt);
        
    end
    
    for i = 1:Nx
        c_exa(i) = c_exact(t,x(i));
    end
    
    error_cn(k) = max(abs(c_exa-c_cn_new'));
    error_im(k) = max(abs(c_exa-c_im_new'));
    
    if k > 1
        order_cn(k) = log(error_cn(k-1)/error_cn(k))/log(2);
        order_im(k) = log(error_im(k-1)/error_im(k))/log(2);
    end
    
    dx_all(k) = dx;
    Nx = Nx*2;
end

fprintf('-----------------------------------------------------------------\n');
fprintf('\t \t Implicit \t\t Crank-Nicolson \n');
fprintf('-----------------------------------------------------------------\n');
fprintf('dx \t\t Error \t\t Order \t Error \t\t Order\n');
fprintf('-----------------------------------------------------------------\n');
for k = 1:num_splits
    fprintf('%g \t %g \t %.2f \t %g \t %.2f \n', dx_all(k), error_im(k), order_im(k), error_cn(k), order_cn(k));
end
fprintf('-----------------------------------------------------------------\n');

