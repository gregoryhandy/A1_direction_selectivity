%%
% The serves as the traveling wave function of an auditory sweep
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
function [b] = spatial_input_fn(x_pos, f_0,b_amp,sweep_speed,curr_t,sigma)

b = b_amp*exp(-(x_pos-log2(f_0*2.^(sweep_speed*curr_t))).^2/sigma^2);

end

