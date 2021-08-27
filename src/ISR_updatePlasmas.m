%%INFO: ISR_main.m���Ӻ���, ��������parameters��plasmas���֡�
%%----------------------------------------------------------------------%%
% Needs: constants.m; 
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - ������̲������м���� struct
% Output:
%   parameters  - ������̲������м���� struct
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/09/10,11
%%----------------------------------------------------------------------%%
%��������[Washy 2020/05/01]
% 1. ȥ������: name, density, temperatur, constants
% 2. �򻯼������, ȥ��if����forѭ��
%%----------------------------------------------------------------------%%

function parameters = ISR_updatePlasmas(parameters)

res  = constants({'kB', 'eps0'});
kB   = res.kB;
eps0 = res.eps0;

omega = parameters.omega;

k = parameters.radar.k;

B = parameters.factors.B;

%% ���Ӳ���

ude = parameters.factors.ude;

me  = parameters.plasmas.me;
qe  = parameters.plasmas.qe;
ne  = parameters.plasmas.ne;
Te  = parameters.plasmas.Te;

% �������ٶ� [m/s]
parameters.plasmas.vTe      = sqrt(2*kB*Te/me);
% ���ӵ�������Ƶ�� [Hz]
parameters.plasmas.wpe      = sqrt(ne*qe^2/eps0/me);
% ���ӵ°ݰ뾶 [m]
parameters.plasmas.he       = sqrt(eps0*kB*Te/ne/qe^2);
% ���ӻ���Ƶ�� [Hz]
parameters.plasmas.Oe       = qe*B/me;
% ���Ӷ�����Ƶ�� [rad*Hz]
parameters.plasmas.omegae   = omega - k*ude;

%% ���Ӳ���

udi = parameters.factors.udi;

mi  = parameters.plasmas.mi;
qi  = parameters.plasmas.qi;
ni  = parameters.plasmas.ni;
Ti  = parameters.plasmas.Ti;

% �������ٶ� [m/s]
parameters.plasmas.vTi      = sqrt(2*kB*Ti./mi);
% ���ӵ�������Ƶ�� [Hz]
parameters.plasmas.wpi      = sqrt(ni.*qi.^2/eps0./mi);
% ���ӵ°ݰ뾶 [m]
parameters.plasmas.hi       = sqrt(eps0*kB.*Ti./ni./qi.^2);
% ���ӻ���Ƶ�� [Hz]
parameters.plasmas.Oi       = qi*B./mi;
% ���Ӷ�����Ƶ�� [rad*Hz]
parameters.plasmas.omegai   = omega - k*udi;

end