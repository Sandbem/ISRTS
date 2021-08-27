%%INFO: ����ISR_main
% date: 2020-05-01 Washy

clear;

ion = {'O+'};
ne = 1e12;
Te = 1000;
Ti = 1000;
percent = 1;

lfre = 1000;
frequency = linspace(0, 2000, lfre);

fradar = 50e6;
theta = 180;

factors = struct;
factors.ud = [0, 0];
factors.nu = [0, 0];
factors.B = [3.5e-5, 87];
factors.mode = [0, 0, 1, 0];
factors.gordmode = 2;
factors.theomode = 1;

[spec,parameters]=ISR_main(ion,ne,Te,Ti,percent,frequency,fradar,theta,factors);

%%
figure;

plot(frequency, spec, 'LineWidth', 1.5)
grid on;