%% 
% Calculates the DSI given the upward and downward responses
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
function [ DSI, area_up,area_down ] = DSI_calc(t,x, delta_r_up,delta_r_down)

%% Subtract out baseline firing rate
delta_r_up = delta_r_up-delta_r_up(1,:);
delta_r_down = delta_r_down-delta_r_down(1,:);

%% Calculate the area under the curve
area_up = zeros(length(x),1);
area_down = zeros(length(x),1);
for i = 1:length(x)
    area_up(i) = trapz(t,delta_r_up(:,i));
    area_down(i) = trapz(t,delta_r_down(:,i));
end

DSI = (area_up-area_down)./(area_up+area_down);

end

