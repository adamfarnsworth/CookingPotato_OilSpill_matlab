clear;
clc;

% domain
xL = -2;
xR = 2;

% range
yB = -2.5;
yT = 2.5;

% diffusion coefficient
D = 0.0015;

% time interval
t_start = 0;
t_final = 1500;

% time-step 
dt = 5;

% source term
f = @(t,x,y) 0;

% initial condition
c_start = @(x,y) 20;

% boundary condition
c_bc = @(t,x,y) min(20+80*t/60,100);

% space discretization
Nx = 80;
Ny = 100;

    dx = (xR-xL)/(Nx-1);
    dy = (yT-yB)/(Ny-1);

    % Create x and y arrays
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
%     x = linspace(xL, xR, Nx);
%     y = linspace(yB, yT, Ny);
%     
%     x_temp = x;
%     y_temp = y;
%    for j = 1:Ny-1
%       x = [x,x_temp];
%    end
%    
%    for j = 1:Nx-1
%       y = [y,y_temp];
%    end


        % Impose initial conditions
        for i=1:Nx
            for j=1:Ny
                c_old(j,i) = c_start(x(i),y(j));
            end
        end
        c_new = c_old;


    t = t_start;

    % Create sparse matrix and allocate memory for right-hand side
    A = sparse(Nx*Ny,Nx*Ny);
   % A = eye(Nx*Ny,Nx*Ny);
    RHS = zeros(Nx*Ny,1);

    %a's
    aC = 1+2*dt*D/dx/dx+2*dt*D/dy/dy;
    aL = -dt*D/dx/dx;
    aR = -dt*D/dx/dx;
    aB = -dt*D/dy/dy;
    aT = -dt*D/dy/dy;

    % Calculate the matrix before the while-loop to save time
    % internal points


    for i = 1:Nx
        for j = 1:Ny
        if ((i == 1) || (i == Nx) || (j == 1) || (j==Ny))
            A((j-1)*Nx+i,j)= 0;
            A((j-1)*Nx+i,(j-1)*Nx+i)= 1;
        else
        A((j-1)*Nx+i,(j-1)*Nx+i) = aC;
        A((j-1)*Nx+i,(j-1)*Nx+i-1) = aL;
        A((j-1)*Nx+i,(j-1)*Nx+i+1) = aR;
        A((j-1)*Nx+i,(j-1)*Nx+i-Nx) = aB;
        A((j-1)*Nx+i,(j-1)*Nx+i+Nx) = aT;
        end
     
        end
    end

    % boundary points
    A(1,1) = 1;
    A(Nx*Ny,Nx*Ny) = 1;
    
% Find size of arrays
n_total = ceil((t_final - t_start)/dt);

% Create preallocated t and potato
t = zeros(1, n_total);
potato = zeros(1, n_total);

% Put initial values into t and potato
t(1) = t_start;
potato(1) = 20;
subplot_num = 1;
% Initialize auxiliary variable to keep track of the number of iterations
n = 1;
figure('rend', 'painters', 'pos', [400 50 600 900])
    while t(n) < t_final
        if((t(n) == 0) ||(t(n) == 200) || (t(n) == 400) || (t(n) == 600))
           subplot(4,2,subplot_num)
           contourf(x,y,c_new, 100, 'LineColor', 'none');
           title(['Temperature of potato at t = ' num2str(t(n))]);
           colorbar;
           caxis([20,100]);
           hold on
           subplot_num = subplot_num +1;
        end
        if t(n) + dt > t_final
            dt = t_final-t(n);
            % need to recalculate the matrix since dt has changed
            % internal points
                aC = 1+2*dt*D/dx/dx+2*dt*D/dy/dy;
                aL = -dt*D/dx/dx;
                aR = -dt*D/dx/dx;
                aB = -dt*D/dy/dy;
                aT = -dt*D/dy/dy;

           for i = 1:Nx
                for j = 1:Ny
                if ((i == 1) || (i == Nx) || (j == 1) || (j==Ny))
                    A((j-1)*Nx+i,j)= 0;
                    A((j-1)*Nx+i,(j-1)*Nx+i)= 1;
                else
                A((j-1)*Nx+i,(j-1)*Nx+i) = aC;
                A((j-1)*Nx+i,(j-1)*Nx+i-1) = aL;
                A((j-1)*Nx+i,(j-1)*Nx+i+1) = aR;
                A((j-1)*Nx+i,(j-1)*Nx+i-Nx) = aB;
                A((j-1)*Nx+i,(j-1)*Nx+i+Nx) = aT;
                end

                end
            end
        end
          %c_old = reshape(c_old.',[],1);
        
         for i = 1:Nx
            for j = 1:Ny
                 if ((i == 1) || (i == Nx) || (j == 1) || (j==Ny))
                     RHS((j-1)*Nx+i) = c_bc(t(n)+dt,x(i),y(j));                      
                 else
                      RHS((j-1)*Nx+i) = c_old(j,i) + dt*f(t(n)+dt,x(i),y(j));
                 end
            end
         end

        c_new = vec2mat(A\RHS,Nx);

        c_old = c_new;
        % Calculate new time
        t(n+1) = t(n) + dt;
        
        % Calculate new potato temp in center
        potato(n+1) = c_new(Ny/2,Nx/2);
        n=n+1;

    end
%                     % Plot 
        subplot(4,2,[5,8]);
        plot(t,potato);
        xlabel('Time');
        ylabel('Temp');
        %A = full(A);



