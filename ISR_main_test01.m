%%INFO: ����ISR_main
% date: 2020-05-01 Washy

clear;

ion = {'O+', 'NO+'};
ne = 5e10;
Te = 1000;
Ti = 1000;
percent = [0.7 0.3];

lfre = 1000;
frequency = linspace(0, 100e3, lfre);

fradar = 440e6;
theta = 180;

factors = struct;
factors.ud = [0, 0];
factors.nu = [0, 5000, 5000];
factors.B = [3.5e-5, 60];
factors.mode = [0, 0, 1, 0];
factors.gordmode = 1;
factors.theomode = 2;

[spec,parameters]=ISR_main(ion,ne,Te,Ti,percent,frequency,fradar,theta,factors);

%%
figure;

plot(frequency, spec, 'LineWidth', 1.5)
grid on;