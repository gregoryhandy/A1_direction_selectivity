%%
% Quick solves for the steady state solution of the non-spatial model
% Used when adjusting parameters assoicated with SOM and PV block
%
% This code was writen by Gregory Handy (2020)
% Please email ghandy@uchicago.edu with any questions
%%
function [ drdt ] = A1_quick_SS_ODE(t, r, param)

curr = param.W*r + param.inhibitory_curr';

drdt = (-r + curr)./param.tau_A';
% thresholding for E population
if r(1)<param.threshold && drdt(1)<0
    drdt(pc)=0;
end

end