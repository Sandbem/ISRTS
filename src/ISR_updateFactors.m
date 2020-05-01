%%INFO: ISR_main.m���Ӻ���, ��������parameters��factors���֡�
%%----------------------------------------------------------------------%%
% Needs: CollisionCoulomb.m
%%----------------------------------------------------------------------%%
% Inputs:
%   parameters  - ������̲������м���� struct
%   factors     - Ư�ơ���ײ���ų���ģʽ struct
% Outputs:
%   parameters  - ������̲������м���� struct
%   mode        - ģʽ����(��: 0; ��: 1)
%   theomode    - �����׼���ģʽ(Mac: 1; KM: 2)
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/09/10,11
%%----------------------------------------------------------------------%%
%��������[Washy 2019/09/24]
%%INFO: ����ײƵ��nu�ӵ�����ֵ�����Ϊ���������룬���е�һ����ֵ��ʾ����
%   ��ײƵ�ʣ��ڶ��������һ����ʾ������ײƵ�ʡ�
% 1. �޸���nu�����ݳ��ȣ�
% 2. �޸���if���mode(2)�еĸ�ֵ��ʽ��
%%----------------------------------------------------------------------%%
%��������[Washy 2019/10/08]
%%INFO: ��ISRKM_analysis_factors.m�ϲ����˸ú��������Ժ���������ṹ����
%   �����¹滮��
% 1. ȥ���������: ud, nu, B, mode
% 2. �����������: factors
% 3. �����������: mode, gordmode
% 4. �޸ı���factors.ud: ���ֵ��ӡ��������룬�޸ĺ�ĸ�ʽΪ[ude, udi]��
%   udiֻ��һ��ֵ����Ĭ�����гɷ־�����ͬ��ֵ��
%%----------------------------------------------------------------------%%
%��������[Washy 2019/10/10]
% 1. �����������: theomode
%%----------------------------------------------------------------------%%
%��������[Washy 2020/05/01]
% 1. ȥ������: ion, ne, ni, Te, Ti
% 2. ȥ�����: gordmode
%%----------------------------------------------------------------------%%

function [parameters, mode, theomode] = ISR_updateFactors(parameters,factors)

leni = length(parameters.plasmas.namei);

% ȱʡֵ
B = [3.5e-5, 0];

% ģʽ����(��: 0; ��: 1)
if isfield(factors,'mode')
    validateattributes(factors.mode,{'double'},{'row','binary'})
    mode = factors.mode;
    
    % Ư���ٶ�
    if isfield(factors,'ud') && mode(1)
        parameters.factors.ude    = factors.ud(1);
        parameters.factors.udi    = factors.ud(2);
    else
        parameters.factors.ude    = 0;
        parameters.factors.udi    = 0;
    end
    
    % ��ײƵ��
    if isfield(factors,'nu') && mode(2)
        parameters.factors.nuen = factors.nu(1);
        parameters.factors.nuin = factors.nu(2:end);
    else
        parameters.factors.nuen = 0;
        parameters.factors.nuin = zeros(1,leni);
    end
    
    % �ų�
    if isfield(factors,'B') && mode(3)
        parameters.factors.B     = factors.B(1);
        parameters.factors.alpha = factors.B(2);
    else
        parameters.factors.B     = B(1);
        parameters.factors.alpha = B(2);
    end
    
    % ������ײƵ��
    if mode(4)
        [parameters.factors.nuii, parameters.factors.nuei, ...
            parameters.factors.nuee] = ...
            CollisionCoulomb(parameters.plasmas.namei, ...
            parameters.plasmas.Ti, parameters.plasmas.Te, ...
            parameters.plasmas.ni, parameters.plasmas.ne);
    else
        parameters.factors.nuii = zeros(1,leni);
        parameters.factors.nuei = 0;
        parameters.factors.nuee = 0;
    end
    
else
    mode = [0,0,0,0];
    
    parameters.factors.ude   = 0;
    parameters.factors.udi   = 0;
    
    parameters.factors.nuen  = 0;
    parameters.factors.nuin  = zeros(1,leni);
    
    parameters.factors.B     = B(1);
    parameters.factors.alpha = B(2);
    
    parameters.factors.nuii  = zeros(1,leni);
    parameters.factors.nuei  = 0;
    parameters.factors.nuee  = 0;
    
end

parameters.factors.mode = mode;

% �����׼��㹫ʽѡ��
if isfield(factors,'theomode')
    validateattributes(factors.theomode,{'double'},{'scalar','integer','>=',1,'<=',2})
    theomode = factors.theomode;
else
    theomode = 1;
end

end