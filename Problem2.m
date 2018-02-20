clear;
clc;


% domain
xL = -1;
xR = 3;

% range
yB = -1.5;
yT = 1.5;

% diffusion coefficient
D = 0.7;

% time interval
t_start = 0;
t_final = 1;

%Velocity
Vx = -0.8;
Vy = -0.4;

% exact solution
c_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% source term  vx -vy
%  f = @(t,x,y) exp(-t)*((2*D-t)*sin(x)*cos(y)+Vx*cos(x)*cos(y)-Vy*sin(y)*sin(x));
 f = @(t,x,y) exp(-t)*(-sin(x)*cos(y) + Vx*cos(x)*cos(y)-Vy*sin(x)*sin(y)+2*D*sin(x)*cos(y));

% initial condition
c_start = @(x,y) c_exact(t_start,x,y);

% boundary condition
c_bc = @(t,x,y) c_exact(t,x,y);

%that g function
g = @(t,x,y) -D*sin(y)*sin(x)*exp(-t)-Vy*sin(x)*cos(y)*exp(-t);

%boundary conditoin bottom and corners taken care of by robbin condition

% space discretization
Nx = 20;
Ny = 15;

num_splits = 4;
for k = 1:num_splits

    dx = (xR-xL)/(Nx-1);
    dy = (yT-yB)/(Ny-1);

    % Create x and y arrays
%     x = linspace(xL, xR, Nx*Ny);
%     y = linspace(yB, yT, Nx*Ny);
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
    
    x_temp = x;
    y_temp = y;
   for j = 1:Ny-1
      x = [x,x_temp];
   end
   
   for j = 1:Nx-1
      y = [y,y_temp];
   end
    % time-step 
    dt = 0.5*dx;

        % Impose initial conditions
        for i=1:Nx
            for j=1:Ny
                c_old(j,i) = c_start(x(i),y(j));
            end
        end
        c_new = c_old;

    c_exa = zeros(1,Nx);


    t = t_start;

    % Create sparse matrix and allocate memory for right-hand side
%     A = zeros(Nx*Ny,Nx*Ny);
    A = sparse(Nx*Ny,Nx*Ny);
    %A = eye(Nx*Ny,Nx*Ny);
    RHS = zeros(Nx*Ny,1);

    %a's
    aC = 1+2*dt*D/dx/dx+2*dt*D/dy/dy;
    aL = -dt*D/dx/dx;
    aR = -dt*D/dx/dx;
    aB = -dt*D/dy/dy;
    aT = -dt*D/dy/dy;
    a_L= -dt*D/dx/dx;
    a_R= -dt*D/dx/dx;
    a_C= aC + 2*Vy*dt/dy;
    a_T= aT + aB;

    % Calculate the matrix before the while-loop to save time
    % internal points


    for i = 1:Nx
        for j = 1:Ny
        if ((i == 1) || (i == Nx) || (j == 1) || (j==Ny))
            A((j-1)*Nx+i,j)= 0;
            A((j-1)*Nx+i,(j-1)*Nx+i)= 1;
            
            if((j == 1) && ((i ~= 1) && (i ~= Nx)) )
                A((j-1)*Nx+i,(j-1)*Nx+i)= a_C;
                A((j-1)*Nx+i,(j-1)*Nx+i + 1)= a_R;
                A((j-1)*Nx+i,(j-1)*Nx+i - 1)= a_L;
                A((j-1)*Nx+i,(j-1)*Nx+i + Nx)= a_T;
            end
                
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

    while t < t_final

        if t + dt > t_final
            dt = t_final-t;
            % need to recalculate the matrix since dt has changed
            % internal points
                   aC = 1+2*dt*D/dx/dx+2*dt*D/dy/dy;
                    aL = -dt*D/dx/dx;
                    aR = -dt*D/dx/dx;
                    aB = -dt*D/dy/dy;
                    aT = -dt*D/dy/dy;
                    a_L= -dt*D/dx/dx;
                    a_R= -dt*D/dx/dx;
                    a_C= aC + 2*Vy*dt/dy;
                    a_T= aT + aB;

           for i = 1:Nx
                for j = 1:Ny
                if ((i == 1) || (i == Nx) || (j == 1) || (j==Ny))
                    A((j-1)*Nx+i,j)= 0;
                    A((j-1)*Nx+i,(j-1)*Nx+i)= 1;
                     
                    if((j == 1) && ((i ~= 1) && (i ~= Nx)) )
                        A((j-1)*Nx+i,(j-1)*Nx+i)= a_C;
                        A((j-1)*Nx+i,(j-1)*Nx+i + 1)= a_R;
                        A((j-1)*Nx+i,(j-1)*Nx+i - 1)= a_L;
                        A((j-1)*Nx+i,(j-1)*Nx+i + Nx)= a_T;
                    end             
                    
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
                     RHS((j-1)*Nx+i) = c_bc(t+dt,x(i),y(j));          
                     
                        if((j == 1) && ((i ~= 1) && (i ~= Nx)) )
                       RHS((j-1)*Nx+i) = c_old(j,i) + dt*f(t+dt,x(i),y(j))...
                                      -Vx*dt*(c_old(j,i+1) - c_old(j,i))/dx...
                                      -Vy*dt*(c_old(j+1,i) - c_old(j,i))/dy ...
                                      -2*dt/dy*g(t+dt,x(i),y(1));
                        end 
                     
                 else
                      RHS((j-1)*Nx+i) = c_old(j,i) + dt*f(t+dt,x(i),y(j))...
                                      -Vx*dt*(c_old(j,i+1) - c_old(j,i))/dx...
                                      -Vy*dt*(c_old(j+1,i) - c_old(j,i))/dy ;
                 end
            end
         end

        c_new = vec2mat(A\RHS,Nx);

        c_old = c_new;
        t = t+dt;

        % internal points
        for i = 1:Nx
            for j=1:Ny
            c_exa(j,i) = c_exact(t,x(i),y(j));
            end
        end
        
        c_new_error = reshape(c_new.',[],1);
        c_exa_error = reshape(c_exa.',[],1);
        error(k) = max(abs(c_exa_error-c_new_error));      
    end
    if k > 1
         order(k) = log(error(k-1)/error(k))/log(2);
    fprintf('Matrix size = %i X %i \nerror = %i \norder of accuracy = %1.5f\n\n', Nx, Ny, error(k), order(k));
    end
    Nx = Nx*2;
    Ny = Ny*2;
end

