clear;
clc;
clf;


oneCount = 0;
fourCount = 0;
sevenCount = 0;

% domain
xL = 0;
xR = 12;

% range
yB = 0;
yT = 3;

% diffusion coefficient
D = 0.2;

% time interval
t_start = 0;
t_final = 10;

%Velocity
Vx = -0.8;
Vy = -0.4;

% % source term  vx -vy
%  f = @(t,x,y) 1/2*( 1-tanh((sqrt( (x-Xs)^2+y^2)-rs )/e) );

% initial condition
c_start = @(x,y) 0;

% boundary condition
c_bc = @(t,x,y) 0;

%that g function
g = @(t,x,y) 0;

%boundary conditoin bottom and corners taken care of by robbin condition

% space discretization
Nx = 160;
Ny = 40;
beach = 4;

num_splits = 3;
for k = 1:num_splits


    dx = (xR-xL)/(Nx-1);
    dy = (yT-yB)/(Ny-1);

    % Create x and y arrays
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
%     x = linspace(xL, xR, Nx);
%     y = linspace(yB, yT, Ny);
    
%     x_temp = x;
%     y_temp = y;
%    for j = 1:Ny-1
%       x = [x,x_temp];
%    end
%    
%    for j = 1:Nx-1
%       y = [y,y_temp];
%    end
    % time-step 
    dt = 0.1;

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

    % Find size of arrays
n_total = ceil((t_final - t_start)/dt);

% Create preallocated t and potato
t = zeros(1, n_total);
oil = zeros(1, n_total);

% Put initial values into t and potato
t(1) = t_start;
oil(1) = 0;
subplot_num = 1;
% Initialize auxiliary variable to keep track of the number of iterations
n = 1;

    while t(n) < t_final


        if t(n) + dt > t_final
            dt = t_final-t(n);
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
                     RHS((j-1)*Nx+i) = c_bc(t(n)+dt,x(i),y(j));          
                     
                        if((j == 1) && ((i ~= 1) && (i ~= Nx)) )
                       RHS((j-1)*Nx+i) = c_old(j,i) + dt*f(t(n)+dt,x(i),y(j))...
                                      -Vx*dt*(c_old(j,i+1) - c_old(j,i))/dx...
                                      -Vy*dt*(c_old(j+1,i) - c_old(j,i))/dy ...
                                      -2*dt/dy*g(t(n)+dt,x(i),y(1));
                        end 
                     
                 else
                      RHS((j-1)*Nx+i) = c_old(j,i) + dt*f(t(n)+dt,x(i),y(j))...
                                      -Vx*dt*(c_old(j,i+1) - c_old(j,i))/dx...
                                      -Vy*dt*(c_old(j+1,i) - c_old(j,i))/dy ;
                 end
            end
         end

        c_new = vec2mat(A\RHS,Nx);

        c_old = c_new;

        
        if ((t(n) > .999) && (t(n) < 1.001) && (oneCount ==0))
         disp('made it to 1 for beach');
         disp(beach);
         contourf(x,y,c_new, 100, 'LineColor', 'none');
         colorbar;
         hold on
         figure; 
         oneCount = oneCount +1;
        end
        
         if ((t(n) > 3.999) && (t(n) < 4.001) && (fourCount ==0))
         disp('made it to 4 for beach');
         disp(beach);
         contourf(x,y,c_new, 100, 'LineColor', 'none');
         colorbar;
         hold on
         figure; 
         fourCount = fourCount +1;
         end
        
         if ((t(n) > 6.999) && (t(n) < 7.001) && (sevenCount ==0))
         disp('made it to 7 for beach');
         disp(beach);
         contourf(x,y,c_new, 100, 'LineColor', 'none');
         colorbar;
         hold on
         figure; 
         sevenCount = sevenCount +1;
        end
        % Calculate new time
        t(n+1) = t(n) + dt;
        
        % Calculate new oil concentration
        oil(n+1) = c_new(1,floor(beach/dx+1));
        disp(x(floor(beach/dx+1)));
        n=n+1;


      
    end
    %disp(beach);
%     A = full(A);

%          contourf(x,y,c_new, 100, 'LineColor', 'none');
%          colorbar;
%          hold on
%          figure;
        plot(t,oil);
        xlabel('Time');
        ylabel('Concentration');
            if (beach < 8)
             figure;
            end
        
        beach = beach +2;
        
end

function hello = f(t,x,y)
%those random variables
Xs = 10;
rs = 0.1;
e = 0.1;
if(t>0.5)
    hello = 0;
else
hello = 1/2*( 1-tanh((sqrt( (x-Xs)^2+y^2)-rs )/e) );
end

end
