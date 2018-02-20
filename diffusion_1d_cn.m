% Crank-Nicolson scheme

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
Nx = 100;
x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);

% time-step from stability
dt = 5*dx*dx/2/D;
%dt = 0.01*dx;

c_old = zeros(1,Nx);
c_new = zeros(1,Nx);
c_exa = zeros(1,Nx);

for i = 1:Nx
    c_old(i) = c_start(x(i));
end

t = t_start;

A = sparse(Nx,Nx);
RHS = zeros(Nx,1);
    
% internal points
for i = 2:Nx-1
    A(i,i) = 1+dt*D/dx/dx;
    A(i,i-1) = -.5*dt*D/dx/dx;
    A(i,i+1) = -.5*dt*D/dx/dx;
end

% boundary points
A(1,1) = 1;
A(Nx,Nx) = 1;
    
while t < t_final
    
    if t + dt > t_final
        dt = t_final-t;
        
        % internal points
        for i = 2:Nx-1
            A(i,i) = 1+dt*D/dx/dx;
            A(i,i-1) = -.5*dt*D/dx/dx;
            A(i,i+1) = -.5*dt*D/dx/dx;
        end
    end
    
    % internal points
    for i = 2:Nx-1        
        RHS(i) = c_old(i) + 0.5*dt*(f(t,x(i))+f(t+dt,x(i))) + .5*dt*D*(c_old(i-1)-2*c_old(i)+c_old(i+1))/dx/dx;
    end
    
    % boundary points
    RHS(1) = cL(t+dt);
    RHS(Nx) = cR(t+dt);
    
    % solve system of equations
    c_new = A\RHS;
    
    c_old = c_new;
    t = t+dt;
    
    for i = 1:Nx
        c_exa(i) = c_exact(t,x(i));
    end
    
    plot(x,c_exa,'LineWidth',2);
    hold on
    plot(x,c_new,'o-','LineWidth',1);
    hold off
    xlabel('x');
    ylabel('c');
    axis([xL xR 0 1]);
    pause(1000*dt);
    
end