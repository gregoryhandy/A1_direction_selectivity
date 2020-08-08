%%
% Recreate Figure 4 using the data produced from A1_run_Fig4and5_sims.m 
% 
% Note: some numerical artifacts in DSI curves might appear due to
% adjustments to the code to optimize speed and sizes of data
%
% This code was writen by Gregory Handy (2020)
% Please email ghandy@uchicago.edu with any questions
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('A1_create_fig4.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Choose the name of the data directory and wave speed

% Default
data_dir = './A1_Sim_Data_Default';
wave_speed = 20*10^-3;
plot_abs_dsi = 1;

% Supplemental
% data_dir = './A1_Sim_Data_SuppFig5';
% wave_speed = 10*10^-3;
% plot_abs_dsi = 0;

% FFWD
% data_dir = './A1_Sim_Data_FFWD';
% wave_speed = 20*10^-3;
% plot_abs_dsi = 0;

temp_file_name = sprintf('/freq_sweep_%.1f.mat',wave_speed*10^3);
load(strcat(data_dir,temp_file_name));

%%  Pick a BF to plot
x_interested = find(2.^param.x > 5,1);

%% Load the incoming excitatory input for the E population at BF specified

[K_up, K_down, recurrent_exc_up, recurrent_exc_down] =...
    recurrent_input(param, x_interested, delta_r_up,...
    delta_r_down);

%% Load upward and downward feeforward input at specific BF

[curr_input_up, curr_input_down, ffwd_total_up, ffwd_total_down] =...
    ffwd_input(param, x_interested);

total_up = ffwd_total_up + recurrent_exc_up;
total_down = ffwd_total_down + recurrent_exc_down;

data_up = [ffwd_total_up,total_up]/ffwd_total_up;
data_down = [ffwd_total_down,total_down]/ffwd_total_down;

%% Create EPSC, total EPSC, etc. figure
x_min = -50;
x_max = param.stim_stop-param.stim_delay+200;
param.color_scheme([2 1],:) = param.color_scheme([1 2],:);

figure(1)
subplot(5,1,1)
hold off
plot(param.tspan-param.stim_delay,curr_input_up','linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan-param.stim_delay,curr_input_down','linewidth',1.5,'color',param.color_scheme(2,:))
title('BF = 5 (kHz)')
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-0.1 11])
ylabel('Stim EPSC')
xticks([])
box off 

subplot(5,1,2)
hold off
plot(param.tspan-param.stim_delay,K_up(1,:)+curr_input_up','linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan-param.stim_delay,K_down(1,:)+curr_input_down','linewidth',1.5,'color',param.color_scheme(2,:))
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-2 20])
ylabel('\Delta Total EPSC')
xticks([])
box off 

subplot(5,1,3)
hold off
plot(param.tspan-param.stim_delay,-K_up(2,:),'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan-param.stim_delay,-K_down(2,:),'linewidth',1.5,'color',param.color_scheme(2,:))
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-3 11])
ylabel('\Delta PV IPSC')
xticks([])
box off 

subplot(5,1,4)
hold off
plot(param.tspan-param.stim_delay,-K_up(3,:),'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan-param.stim_delay,-K_down(3,:),'linewidth',1.5,'color',param.color_scheme(2,:))
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-0.5 4])
ylabel('\Delta SOM IPSC')
xticks([])
box off 

% subplot(5,1,5)
% hold off
% plot(param.tspan-param.stim_delay,-K_up(2,:)-K_up(3,:),'linewidth',1.5,'color',param.color_scheme(1,:))
% hold on
% plot(param.tspan-param.stim_delay,-K_down(2,:)-K_down(3,:),'linewidth',1.5,'color',param.color_scheme(2,:))
% set(gca,'fontsize',16)
% xlim([x_min x_max])
% ylim([-0.5 11])
% ylabel('\Delta Total IPSC')
% xticks([])
% box off 

subplot(5,1,5)
hold off
plot(param.tspan-param.stim_delay,(delta_r_up(:,x_interested)-delta_r_up(1,x_interested))+1,'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan-param.stim_delay,(delta_r_down(:,x_interested)-delta_r_down(1,x_interested))+1,'linewidth',1.5,'color',param.color_scheme(2,:))

ylim([-1 13])
xlim([x_min x_max])
xticks([0:100:x_max])
set(gca,'fontsize',16)
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
box off 
legend('Upward Sweep','Downward Sweep')
legend box off

%% Create the EPSC amplification curve

x_temp = [1, 2];
figure(2)
hold off
plot(x_temp,data_up,'-','markersize',16,'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(x_temp,data_down,'-','markersize',16,'linewidth',1.5,'color',param.color_scheme(2,:))

set(gca,'fontsize',16)
xlim([0.8 2.2])
xticks([1 2])
xticklabels({'Feedforward','Total'})
ylabel('Normalized EPSC')
legend('Upward (preferred)','Downward (non-pref)')
legend box off
box off

%% Create the DSI curve over the tonotopic space

figure(3)
hold off
plot(param.x,DSI_values,'k','linewidth',1.5)
hold on
plot(param.x,param.x*0,'k--','linewidth',1)
plot(param.x(x_interested),DSI_values(x_interested),'.','markersize',20,'color',[0.8 0 0])
ylim([-1 1])
xlim([2.27 5.73])
set(gca,'fontsize',16)
xticks([2 3 4 5 6])
xticklabels({'4','8','16','32','64'})
ylabel('DSI')
xlabel('BF (kHz)')
box off

%% Now, lets load the data for the other wave speeds and create the abs(DSI) panel

if plot_abs_dsi == 1
    sweep_speeds = [2.5 5 10 20 40 80];
    x_index_1 = find(param.x>2.27,1);
    x_index_2 = find(param.x>5.73,1);

    abs_DSI = zeros(length(sweep_speeds),1);
    for ii = 1:length(sweep_speeds)
        temp_file_name = sprintf('/freq_sweep_%.1f.mat',sweep_speeds(ii));
        load(strcat(data_dir,temp_file_name))
        
        abs_DSI(ii) = trapz(param.x(x_index_1:x_index_2),abs(DSI_values(x_index_1:x_index_2)));
    end
    
    figure(4);
    plot([1:length(sweep_speeds)],abs_DSI,'k-','linewidth',1.5)
    set(gca,'fontsize',16)
    xticks([1:length(sweep_speeds)])
    xticklabels({'2.5','5','10','20','40','80'})
    xlabel('Wave Speed (oct/sec)')
    ylabel('Absolute DSI (AUC)')
    ylim([0 1])
end
box off
%%


%%

t_start = find(param.tspan>450,1);
t_end = find(param.tspan>1000,1);

opt_colors =[0.4940    0.1840    0.5560; 0.4660    0.6740    0.1880];

figure(5);
hold off
plot(param.tspan(t_start:t_end)-param.stim_delay,K_down(1,t_start:t_end),'linewidth',1.5,'color',opt_colors(1,:))
hold on
plot(param.tspan(t_start:t_end)-param.stim_delay,-K_down(2,t_start:t_end)-K_down(3,t_start:t_end),'linewidth',1.5,'color',opt_colors(2,:))
% plot(param.tspan(t_start:t_end),param.tspan(t_start:t_end)*0,'k--')
set(gca,'fontsize',16)
legend('EPSC','IPSC')
legend box off
xlim([param.tspan(t_start)-param.stim_delay param.tspan(t_end)-param.stim_delay])
ylabel('Total EPSC or IPSC')
xlabel('Time (ms)')
title('BF = 5 (Hz)')
box off


t_start = find(param.tspan>450,1);
t_end = find(param.tspan>750,1);
axes('Position',[.8 .75 .1 .1])
box on 
plot(param.tspan(t_start:t_end)-param.stim_delay,K_down(1,t_start:t_end),'linewidth',1.5,'color',opt_colors(1,:))
hold on
plot(param.tspan(t_start:t_end)-param.stim_delay,-K_down(2,t_start:t_end)-K_down(3,t_start:t_end),'linewidth',1.5,'color',opt_colors(2,:))
plot(param.tspan(t_start:t_end)-param.stim_delay,param.tspan(t_start:t_end)*0,'k--')
ylim([-1.2 .566]) 
xlim([param.tspan(t_start)-param.stim_delay param.tspan(t_end)-param.stim_delay])
box off