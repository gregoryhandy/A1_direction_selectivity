%%
% Recreate Figure 4 using the data produced from A1_run_all_sims.m 
% 
% Note: some numerical artifacts in DSI curves might appear due to
% adjustments to the code to optimize speed and sizes of data
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('create_fig4.m')); 
addpath(genpath(folder));
rmpath(folder)

%%
load('./A1_Sim_Data/freq_sweep_20_default.mat');

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
data_down = [ffwd_total_down,total_down]/ffwd_total_up;

%%
x_min = 450;
x_max = 850;

figure(1)
subplot(6,1,1)
hold off
plot(param.tspan,curr_input_up','linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan,curr_input_down','linewidth',1.5,'color',param.color_scheme(2,:))
title('BF = 5 (kHz)')
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-0.1 11])
ylabel('Input EPSC')
legend('Upward Sweep','Downward Sweep')

subplot(6,1,2)
hold off
plot(param.tspan,K_up(1,:)+curr_input_up','linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan,K_down(1,:)+curr_input_down','linewidth',1.5,'color',param.color_scheme(2,:))
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-2 20])
ylabel('Total EPSC')

subplot(6,1,3)
hold off
plot(param.tspan,-K_up(2,:),'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan,-K_down(2,:),'linewidth',1.5,'color',param.color_scheme(2,:))
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-3 11])
ylabel('PV IPSC')

subplot(6,1,4)
hold off
plot(param.tspan,-K_up(3,:),'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan,-K_down(3,:),'linewidth',1.5,'color',param.color_scheme(2,:))
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-0.5 4])
ylabel('SOM IPSC')

subplot(6,1,5)
hold off
plot(param.tspan,-K_up(2,:)-K_up(3,:),'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan,-K_down(2,:)-K_down(3,:),'linewidth',1.5,'color',param.color_scheme(2,:))
set(gca,'fontsize',16)
xlim([x_min x_max])
ylim([-0.5 11])
ylabel('Total IPSC')

subplot(6,1,6)
hold off
plot(param.tspan,(delta_r_up(:,x_interested)-delta_r_up(1,x_interested))+1,'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(param.tspan,(delta_r_down(:,x_interested)-delta_r_down(1,x_interested))+1,'linewidth',1.5,'color',param.color_scheme(2,:))

ylim([-1 13])
xlim([x_min x_max])
set(gca,'fontsize',16)
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')

%%

x_temp = [1, 2];
figure(2)
hold off
plot(x_temp,data_up,'*-','markersize',16,'linewidth',1.5)
hold on
plot(x_temp,data_down,'.-','markersize',16,'linewidth',1.5)

set(gca,'fontsize',16)
xlim([0.8 2.2])
xticks([1 2])
xticklabels({'Ffwd','Total'})
legend('Upward','Downward')

%%

figure(3)
plot(param.x,DSI_values,'linewidth',1.5)
hold on
plot(param.x(x_interested),DSI_values(x_interested),'.','markersize',20,'color',[0.8 0 0])
ylim([-1 1])
xlim([2.27 5.73])
set(gca,'fontsize',16)
xticks([2 3 4 5 6])
xticklabels({'4','8','16','32','64'})
ylabel('DSI')
xlabel('BF')

%% Now, lets load the data for the other wave speeds

load('./A1_Sim_Data/freq_sweep_2.5_Default.mat')
DSI_values_2p5 = DSI_values;

load('./A1_Sim_Data/freq_sweep_5_Default.mat')
DSI_values_5 = DSI_values;

load('./A1_Sim_Data/freq_sweep_10_Default.mat')
DSI_values_10 = DSI_values;

load('./A1_Sim_Data/freq_sweep_20_Default.mat')
DSI_values_20 = DSI_values;

load('./A1_Sim_Data/freq_sweep_40_Default.mat')
DSI_values_40 = DSI_values;

load('./A1_Sim_Data/freq_sweep_80_Default.mat')
DSI_values_80 = DSI_values;

load('./A1_Sim_Data/freq_sweep_160_Default.mat')
DSI_values_160 = DSI_values;

x_index_1 = find(param.x>2.27,1);
x_index_2 = find(param.x>5.73,1);
abs_DSI(1) = trapz(param.x(x_index_1:x_index_2),abs(DSI_values_2p5(x_index_1:x_index_2)));
abs_DSI(2) = trapz(param.x(x_index_1:x_index_2),abs(DSI_values_5(x_index_1:x_index_2)));
abs_DSI(3) =trapz(param.x(x_index_1:x_index_2),abs(DSI_values_10(x_index_1:x_index_2)));
abs_DSI(4) =trapz(param.x(x_index_1:x_index_2),abs(DSI_values_20(x_index_1:x_index_2)));
abs_DSI(5) =trapz(param.x(x_index_1:x_index_2),abs(DSI_values_40(x_index_1:x_index_2)));
abs_DSI(6) =trapz(param.x(x_index_1:x_index_2),abs(DSI_values_80(x_index_1:x_index_2)));
abs_DSI(7) =trapz(param.x(x_index_1:x_index_2),abs(DSI_values_160(x_index_1:x_index_2)));

figure(4);
plot([1:7],abs_DSI,'.-','markersize',16)
set(gca,'fontsize',16)
xticks([1:7])
xticklabels({'2.5','5','10','20','40','80','160'})
xlabel('Wave Speed (oct/sec)')
ylabel('Absolute DSI (AUC)')

