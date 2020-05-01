%%INFO: ≤‚ ‘ISR_main
% date: 2020-05-01 Washy

clear;

addpath('./src')

ion = {'O+', 'NO+'};
ne = 5e10;
Ti = 1000;
Tr = 1;
percent = [0.7 0.3];

lfre = 1000;
frequency = linspace(-100e3, 100e3, lfre);

fradar = 440e6;
theta = 180;

factors = struct;
factors.ud = [0, 0];
factors.nu = [0, 5000, 5000];
factors.B = [3.5e-5, 60];
factors.mode = [0, 1, 0, 0];
factors.theomode = 1;

[spec,parameters]=ISR_main(ion,ne,Ti,Tr,percent,frequency,fradar,theta,factors);

%%
figure;

plot(frequency, spec, 'LineWidth', 1.5)
grid on;