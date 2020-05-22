%%
% Spatial rate model for an arbitrary number of populations
%
% Tonotopic map: not a periodic boundary
%
% Tested for three populations: E, PV, and SOM, but generalizable for more
% 
% drdt is a (Nx*Npop) x 1 vector
% 1:Nx E where the index denotes x location
% (Nx+1):2Nx P
% (2Nx+1):3Nx S
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
function [ drdt ] = A1_spatial_sweep_6_ODE(t, r, x, dx, Nx, Npop, b_amp, ...
    W, tau_A, sigma, lambda, scaling_factor,sweep_speed,f_0,stim_stop,...
    stim_delay, tau_S, inhibitory_curr,threshold)

drdt = zeros(Nx*Npop*2,1);

% update the filters/auxiliary variables
for pc = 1:Npop
    if tau_S(pc)~=0
        drdt((1+Nx*(2*pc-1)):(2*Nx*pc)) = ...
            (r((1+2*Nx*(pc-1)):(Nx*(2*pc-1)))-r((1+Nx*(2*pc-1)):(2*Nx*pc)))/tau_S(pc);
    else
        % take the other filters to be instantaneous
        drdt((1+Nx*(2*pc-1)):(2*Nx*pc)) = 0;
        r((1+Nx*(2*pc-1)):(2*Nx*pc)) =  r((1+2*Nx*(pc-1)):(Nx*(2*pc-1)));
    end
end

% Matrix of connection strenths (updated at each x position)
K = zeros(Npop,Npop);

% loop through the locations
for i = 1:Nx
    % loop through the populations
    for pc = 1:Npop
        % loop through the in-coming connections 
        for j = 1:Npop
            % check to see if there are spatial projections in-coming
            if lambda(pc,j)~=0
                tempY = W(pc,j)*exp(-abs(x(i)-x).^2/lambda(pc,j)^2)/scaling_factor(pc,j,i).*r(1+Nx*(2*j-1):2*Nx*(j))';
                K(pc,j) = trapz(tempY)*dx;
            else
                K(pc,j) = W(pc,j)*r(i+Nx*(2*j-1));
            end
            
        end
        
        stim_input = 0;
        % delay of stim_delay ms before input
        if t > stim_delay && t < stim_stop && pc~=3
            stim_input = spatial_input_fn(x(i),f_0,b_amp(pc),sweep_speed,...
                t-stim_delay,sigma(pc));
        end
        total_curr = (sum(K(pc,:))+stim_input+inhibitory_curr(pc));
        

        temp_curr = 1/tau_A(pc)*(-r(i+2*Nx*(pc-1))+total_curr);
        
        % rate equation
        drdt(i+2*Nx*(pc-1)) =temp_curr;   
        
        % thresholding for E population
        if pc == 1 && r(i+Nx*(pc-1)) < threshold && temp_curr < 0
            drdt(i+Nx*(pc-1))=0;
        end
    end
end

end

