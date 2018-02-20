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
dt = 0.9*dx*dx/2/D;

% uncomment to get an unstable solution 
% dt = 5*dx*dx/2/D;

c_old = zeros(1,Nx);
c_new = zeros(1,Nx);
c_exa = zeros(1,Nx);

for i = 1:Nx
    c_old(i) = c_start(x(i));
end

t = t_start;

while t < t_final
    
    if t + dt > t_final
        dt = t_final-t;
    end
    
    for i = 2:Nx-1
        c_new(i) = c_old(i) + dt*D*(c_old(i-1)-2*c_old(i)+c_old(i+1))/dx/dx + dt*f(t,x(i));
    end
    
    c_new(1)  = cL(t);
    c_new(Nx) = cR(t);
    
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
    pause(100*dt);
    
end