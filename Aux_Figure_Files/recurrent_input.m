%%
% Calculate the recurrent input to the circuit
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
function [ K_up, K_down, recurrent_exc_up, recurrent_exc_down ] = ...
    recurrent_input(param, x_interested, delta_r_up, delta_r_down)

pc = 1;

K_up = zeros(param.Npop,length(param.tspan));
K_down = zeros(param.Npop,length(param.tspan));

% loop through time
for tt = 1:length(param.tspan)
    
    % loop through incoming projections from each population
    for j = 1:param.Npop
        % E and PV don't have a temporal filter (first case), while SOM do
        if j < 3
            tempY = param.W(pc,j)*exp(-abs(param.x(x_interested)-param.x).^2/param.lambda(pc,j)^2)/param.scaling_factor(pc,j,x_interested).*delta_r_up(tt,1+2*param.Nx*(j-1):param.Nx*(2*j-1));
            K_up(j,tt) = trapz(tempY)*param.dx;
            
            tempY = param.W(pc,j)*exp(-abs(param.x(x_interested)-param.x).^2/param.lambda(pc,j)^2)/param.scaling_factor(pc,j,x_interested).*delta_r_down(tt,1+2*param.Nx*(j-1):param.Nx*(2*j-1));
            K_down(j,tt) = trapz(tempY)*param.dx;
        else
            tempY = param.W(pc,j)*exp(-abs(param.x(x_interested)-param.x).^2/param.lambda(pc,j)^2)/param.scaling_factor(pc,j,x_interested).*delta_r_up(tt,(1+param.Nx*(2*j-1)):(2*param.Nx*j));
            K_up(j,tt) = trapz(tempY)*param.dx;
            
            tempY = param.W(pc,j)*exp(-abs(param.x(x_interested)-param.x).^2/param.lambda(pc,j)^2)/param.scaling_factor(pc,j,x_interested).*delta_r_down(tt,(1+param.Nx*(2*j-1)):(2*param.Nx*j));
            K_down(j,tt) = trapz(tempY)*param.dx;
        end
    end
end

recurrent_exc_up = trapz(param.tspan,K_up(1,:));
recurrent_exc_down = trapz(param.tspan,K_down(1,:));

end

