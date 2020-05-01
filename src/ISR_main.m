%%INFO: ������Ϊ�����׼��㺯��������
%%----------------------------------------------------------------------%%
% Needs: ISR_init.m; 
%   ISR_updateFactors.m; ISR_updatePlasmas; ISR_updateDimensionless.m; 
%   ISR_gordeyev.m; ISR_spectrumKM.m; ISR_spectrumMac.m;
%%----------------------------------------------------------------------%%
% Inputs:
%   ion         - ���ӳɷ� cell
%   ne          - �������ܶ� [m^-3]
%   Ti          - �����¶� [K]
%   Tr          - ���������¶ȱ�
%   percent     - ���ӳɷֱ���
%   frequency   - Ƶ�ʷ�Χ [Hz]
%   fradar      - �״�Ƶ�� [Hz]
%   theta       - ɢ��� [��]
%   factors     - Ư�ơ���ײ���ų������ء�ģʽ struct
% Output:
%   spec        - ������ [s]
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/09,30
%%----------------------------------------------------------------------%%

function [spec,parameters]=ISR_main(ion,ne,Ti,Tr,percent,frequency,fradar,theta,factors)
%% ������
if nargin == 8
    factors = struct;
elseif nargin == 9
    validateattributes(factors,{'struct'},{'nonempty'})
elseif nargin < 8
    error('�����������Ϊ8����')
else
    error('�����������Ϊ9����')
end

validateattributes(ion,{'cell'},{'row'})
validateattributes(ne,{'double'},{'real','>',0})
validateattributes(Ti,{'double'},{'real','>',0})
validateattributes(Tr,{'double'},{'real','>',0})
validateattributes(percent,{'double'},{'row','>=',0,'<=',1})
validateattributes(frequency,{'double'},{'row','real'})
validateattributes(fradar,{'double'},{'real','>',0})

%% ��ʼ��parameters
parameters = ISR_init(ion,ne,Ti,Tr,percent,frequency,fradar,theta);

%% ��������
[parameters, mode, theomode] = ISR_updateFactors(parameters,factors);

%% �����������
parameters = ISR_updatePlasmas(parameters);

%% �����ٲ���
parameters = ISR_updateDimensionless(parameters);

%% Gordeyev����
[Je,Ji] = ISR_gordeyev(parameters, mode, theomode);

%% ���ּ��������������
k   = parameters.radar.k;
ni  = parameters.plasmas.ni;
vTe = parameters.plasmas.vTe;
vTi = parameters.plasmas.vTi;
he  = parameters.plasmas.he;

%% ������
if theomode == 1 % K&M
    
    thetae = parameters.dimensionless.thetae;
    thetai = parameters.dimensionless.thetai;
    
    hi   = parameters.plasmas.hi;
    
    spec = ISR_spectrumKM(Je,Ji,thetae,thetai,k,ne,ni,vTe,vTi,he,hi);
 
elseif theomode == 2 % Mac
    
    qi     = parameters.plasmas.qi;
    qe     = parameters.plasmas.qe;
    omegae = parameters.plasmas.omegae;
    omegai = parameters.plasmas.omegai;
    
    psien  = parameters.dimensionless.psien;
    psiin  = parameters.dimensionless.psiin;
    
    eta = ni.*qi.^2/(ne*qe^2);
    mu  = eta*Tr;

    spec = ISR_spectrumMac(Je,Ji,omegae,omegai,k,ne,he,eta,mu,vTe,vTi,psien,psiin);
    
end

end