%%
% The serves as the traveling wave function of an auditory sweep
%
% This code was writen by Gregory Handy (2020)
% Please email ghandy@uchicago.edu with any questions
%%
function [b] = spatial_input_fn(x_pos, f_0,b_amp,sweep_speed,curr_t,sigma,...
    stim_stop,stim_delay,ffwd_factor)

amp_slope = 1./(stim_stop - stim_delay);
b = (b_amp-(b_amp - b_amp*ffwd_factor)*amp_slope*curr_t)*exp(-(x_pos-log2(f_0*2.^(sweep_speed*curr_t))).^2/sigma^2);


%% Cleaner code
% b_amp_max = b_amp
% b_amp_min = b_amp*ffwd_factor
% (b_amp_min-b_amp_max)
% (stim_stop - stim_delay)
% amp_slope = (b_amp_min-b_amp_max)/(stim_stop - stim_delay)
% b = (b_amp+amp_slope*curr_t)*exp(-(x_pos-log2(f_0*2.^(sweep_speed*curr_t))).^2/sigma^2);

end

