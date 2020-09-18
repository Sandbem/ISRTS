%%INFO: ISR_main.m���Ӻ���, ������ʼ��parameters�Լ�����radar���ֵļ��㡣
%%----------------------------------------------------------------------%%
% Needs: constants.m; analysisIon.m; 
%%----------------------------------------------------------------------%%
% Inputs:
%   frequency   - Ƶ�ʷ�Χ [Hz]
%   fradar      - �״�Ƶ�� [Hz]
%   theta       - ɢ��� [��]
% Output:
%   parameters  - ������̲������м���� struct
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/09/10,11
%%----------------------------------------------------------------------%%
%%��������[Washy 2020/01/09]
% 1. ���ӱ���: parameters.dimensionless.He, parameters.dimensionless.Hi
%%----------------------------------------------------------------------%%
%%��������[Washy 2020/05/01]
% 1. ��������: ion, ne, Ti, Tr, percent
%%----------------------------------------------------------------------%%

function parameters = ISR_init(ion,ne,Ti,Tr,percent,frequency,fradar,theta_s)

res = constants({'c','me','e'});

parameters = struct;

% parameters.frequency = frequency;           % ������Ƶ�� [Hz]
parameters.omega = 2*pi*frequency;          % ������Ƶ�� [rad*Hz]

%% �״����
parameters.radar = struct;

parameters.radar.frequency = fradar;        % �״�Ƶ�� [Hz]
parameters.radar.theta     = theta_s;       % ɢ��� [��]
parameters.radar.k         = 4*pi*fradar/res.c*sin(theta_s*pi/360); % ɢ���ʸ [1/m]

%% ��������
parameters.factors = struct;

parameters.factors.ude    = [];             % ����Ư���ٶ�(���߷���)
parameters.factors.udi    = [];             % ����Ư���ٶ�(���߷���)

parameters.factors.nuen   = [];             % ����-������ײƵ�� [Hz]
parameters.factors.nuin   = [];             % ����-������ײƵ�� [Hz]

parameters.factors.B      = [];             % �Ÿ�Ӧǿ�� [T]
parameters.factors.alpha  = [];             % ��ų��н� [��]

parameters.factors.nuii   = [];             % ����-���ӿ�����ײƵ�� [Hz]
parameters.factors.nuei   = [];             % ����-���ӿ�����ײƵ�� [Hz]
parameters.factors.nuee   = [];             % ����-���ӿ�����ײƵ�� [Hz]

% parameters.factors.mode   = [];             % ���ƿ���: 0 ��; 1 ��
parameters.factors.gordmode = [];           % gordeyev����: 1: Sommerfeld; 2: Bessel
parameters.factors.theomode = [];           % �����׼��㹫ʽ: 1: Kudeki & Milla; 2: Farley et al.

%% �����������
parameters.plasmas = struct;

% ���Ӳ���
parameters.plasmas.namee  = {'E-'};         % ������������ cell
parameters.plasmas.me     = res.me;         % �������� [kg]
parameters.plasmas.qe     = res.e;          % ���ӵ���� [C]

parameters.plasmas.ne     = ne;             % �������ܶ� [m^-3]
parameters.plasmas.Te     = Ti*Tr;          % �����¶� [K]

parameters.plasmas.vTe    = [];             % �������ٶ� [m/s]
parameters.plasmas.wpe    = [];             % ���ӵ�������Ƶ�� [Hz]
parameters.plasmas.he     = [];             % ���ӵ°ݰ뾶 [m]
parameters.plasmas.Oe     = [];             % ���ӻ���Ƶ�� [Hz]

parameters.plasmas.omegae = [];             % ���Ӷ�����Ƶ�� [rad*Hz]

% ���Ӳ���
[mi, qi] = analysisIon(ion);
parameters.plasmas.namei  = ion;            % ������������ cell
parameters.plasmas.mi     = mi';            % �������� [kg]
parameters.plasmas.qi     = qi';            % ���ӵ���� [C]

parameters.plasmas.ni     = ne*percent;     % �������ܶ� [m^-3]
parameters.plasmas.Ti     = Ti;             % �����¶� [K]

parameters.plasmas.vTi    = [];             % �������ٶ� [m/s]
parameters.plasmas.wpi    = [];             % ���ӵ�������Ƶ�� [Hz]
parameters.plasmas.hi     = [];             % ���ӵ°ݰ뾶 [m]
parameters.plasmas.Oi     = [];             % ���ӻ���Ƶ�� [Hz]

parameters.plasmas.omegai = [];             % ���Ӷ�����Ƶ�� [rad*Hz]

%% �����ٲ���
parameters.dimensionless = struct;

parameters.dimensionless.thetae     = [];   % �����������Ա���
parameters.dimensionless.thetai     = [];   % �����������Ա���

parameters.dimensionless.psien      = [];   % ����-������������ײƵ��
parameters.dimensionless.psiin      = [];   % ����-������������ײƵ��

parameters.dimensionless.phie       = [];   % ���������ٻ���Ƶ��
parameters.dimensionless.phii       = [];   % ���������ٻ���Ƶ��

parameters.dimensionless.psiec_par  = [];   % ���������ٿ�����ײƵ��(ƽ��)
parameters.dimensionless.psiec_perp = [];   % ���������ٿ�����ײƵ��(��ֱ)
parameters.dimensionless.psiic      = [];   % ���������ٿ�����ײƵ��

% parameters.dimensionless.He         = [];   % ����������
% parameters.dimensionless.Hi         = [];   % ����������

%% �����ײ���
parameters.spectrum = struct;

parameters.spectrum.Je    = [];             % ����Gordeyev����
parameters.spectrum.Ji    = [];             % ����Gordeyev����

parameters.spectrum.spece = [];             % ����ɢ����
parameters.spectrum.speci = [];             % ����ɢ����

end