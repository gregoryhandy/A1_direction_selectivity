%%
% Calculate the feedword input to the circuit
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
function [ curr_input_up, curr_input_down, ffwd_total_up, ffwd_total_down]...
    = ffwd_input( param, x_interested )

curr_input_up = zeros(length(param.tspan),1);
curr_input_down = zeros(length(param.tspan),1);

for j = 1:length(param.tspan)
    if param.tspan(j)>param.stim_delay && param.tspan(j)< param.stim_stop
        curr_input_up(j) = spatial_input_fn(param.x(x_interested),param.f_0_up,...
            param.b_amp(1),-param.sweep_speed_down,param.tspan(j)-param.stim_delay,param.sigma(1));
        
        curr_input_down(j) = spatial_input_fn(param.x(x_interested),param.f_0_down,...
            param.b_amp(1),param.sweep_speed_down,param.tspan(j)-param.stim_delay,param.sigma(1));
    else
        curr_input_up(j) = 0;
        curr_input_down(j) = 0;
    end
end

ffwd_total_up = trapz(param.tspan,curr_input_up');
ffwd_total_down = trapz(param.tspan,curr_input_down');

end

