%%INFO: �������������������ײƵ�ʡ�
%%----------------------------------------------------------------------%%
% Inputs:
%   densties    - ���ܶ� [ne, ni] [m^-3]
%   temperatures- �¶� [Te, Ti] [K]
%   ions        - ������Ϣ cell
%   mode        - ����ģʽ, Ĭ��Ϊ 'e_i'
% Outputs:
%   varargout   - ���, ���� mode �Զ�����
%       �������: ���ֻ��һ������ nue ���� [nuee, nuei]
%       ��������: ���ֻ��һ������ nuii
%       ȫ������: ��������������� [nue, nuii]
%       nuii    - ion-ion������ײƵ�� [Hz]
%       nuei    - ele-ion������ײƵ�� [Hz]
%       nuee    - ele-ele������ײƵ�� [Hz]
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Rober Schunk and Andrew Nagy (2009). Ionospheres: Physics, Plasma
%   Physics, and Chemistry Second edition. [Book]
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/01/02
%%----------------------------------------------------------------------%%
%��������[2019/08/20]
% 1. �޸�ni���뵥λΪm^-3
% 2. �޸�Ti�����벻��ion��Ԫ�ظ����ı�
%%----------------------------------------------------------------------%%
%��������[2019/12/12]
% 1. �޸ĺ�����Ϊ: CollisionCoulomb (ԭ������Ϊ: Coulomb_C)
%%----------------------------------------------------------------------%%
%��������[2020/04/30]
% 1. ȡ����csv�ļ���ȡ, ���IonԪ�����Bst����
%%----------------------------------------------------------------------%%
%��������[2020/08/03][��Ҫ]
%INFO: �Ժ�����������Ż�, ���ڿ���ֻ�������[nuee, nuei]������[nuii]����
%   ��ײƵ��, ��Ϊ��Ҫ�������ڿ���ͬʱ�������ܶȡ��¶ȶ�Ӧ����ײƵ�ʡ�
% 1. �޸�����Ϊ: densties, temperatures, ions, mode 
%   (ԭΪ: ion, Ti, Te, ni, ne)
% 2. �޸����Ϊ: varargout ���� mode �Զ�����
%   (ԭΪ: nuii, nuei, nuee)
%%----------------------------------------------------------------------%%
%��������[2020/10/15]
%INFO: �޸���nuii�����������⡣��123��nust��������У���Ҫ���ȶ�ni(:,ia)
%   ����reshape���ټ��㡣
%%----------------------------------------------------------------------%%
%����ʾ��[2020/08/03]
% ne = 1e10;
% ni = [0.5 0.2 0.3]*ne;
% Ti = 1000;
% Te = 1000;
% denstiy = [ne, ni; 0.1*ne, 0.1*ni; 0.2*ne, 0.2*ni];
% temperature = [Te, Ti; 0.95*Te, 0.95*Ti; 0.9*Te, 0.95*Ti];
% ions = {'O+','H+','NO+'};
% mode = 'e_i';
%%----------------------------------------------------------------------%%

function varargout = CollisionCoulomb(densties, temperatures, ions, mode)

% ������������ж�
if nargin == 2
    ions = {};
    mode = 'e';
elseif nargin == 3
    mode = 'e_i';
elseif nargin < 2 || nargin > 4
    error('��������2������������4��������')
end

mode_e = {'e', 'ele', 'electron', 'E', 'Ele', 'Electron'};
mode_i = {'i', 'ion', 'ions', 'I', 'Ion', 'Ions'};
mode_e_i = {'e_i', 'ele_ion', 'ele_ions', 'E_I', 'Ele_Ion', 'Ele_Ions'};

% ���������ײƵ��
if ismember(mode, [mode_e, mode_e_i])
    % ��ȡ���Ӳ���
    ne = densties(:,1)*1e-6;
    Te = temperatures(:,1);
    
    % ele_ele��ײƵ��
    nuee = 54.5/sqrt(2)*ne./Te.^1.5; %[1] eq(4.145)
    
    % ele_ion��ײƵ��
    % nuei = 54.5*ni*Z'.^2/power(Te, 1.5); %[1] eq(4.144)
    nuei = 54.5*ne./Te.^1.5;

    % ���
    varargout{1} = [nuee, nuei];
end

% ����������ײƵ��
if ismember(mode, [mode_i, mode_e_i])
    % ��ȡ���Ӳ���
    ni = densties(:, size(temperatures,2):end)*1e-6;
    Ti = temperatures(:,end);
    
    % �ɼ�����������
    defaultIons = {'H+','He+','C+','N+','O+','CO+','N2+','NO+','O2+','CO2+'};
    % ��ײƵ��ϵ��: �� - s; �� - t
    Bst = [0.900, 1.140, 1.220, 1.230, 1.230, 1.250, 1.250, 1.250, 1.250, 1.260;
        0.280, 0.450, 0.550, 0.560, 0.570, 0.590, 0.590, 0.600, 0.600, 0.610;
        0.102, 0.180, 0.260, 0.270, 0.280, 0.310, 0.310, 0.310, 0.310, 0.320;
        0.088, 0.160, 0.230, 0.240, 0.250, 0.280, 0.280, 0.280, 0.280, 0.300;
        0.077, 0.140, 0.210, 0.220, 0.220, 0.250, 0.250, 0.260, 0.260, 0.270;
        0.045, 0.085, 0.130, 0.140, 0.150, 0.170, 0.170, 0.170, 0.180, 0.190;
        0.045, 0.085, 0.130, 0.140, 0.150, 0.170, 0.170, 0.170, 0.180, 0.190;
        0.042, 0.080, 0.120, 0.130, 0.140, 0.160, 0.160, 0.160, 0.170, 0.180;
        0.039, 0.075, 0.120, 0.120, 0.130, 0.150, 0.150, 0.160, 0.160, 0.170;
        0.029, 0.055, 0.090, 0.090, 0.100, 0.120, 0.120, 0.120, 0.120, 0.140];
    
    % ion_ion��ײƵ��
    [~, ia, ib] = intersect(ions, defaultIons, 'stable');
    lenTi = length(Ti);
    lenia = length(ia);
    lenib = length(ib);
    nuii = zeros(lenTi, length(ions));
%     for iTi = 1:lenTi
%         nust = Bst(ib, ib).*repmat(ni(iTi, ia), length(ia), 1)/Ti(iTi)^1.5; %[1] eq(4.143)
%         nuii(iTi, ia) = sum(nust, 2);
%     end
    nust = repmat(reshape(Bst(ib, ib),1,lenib,lenib), lenTi, 1, 1)...
        .*repmat(reshape(ni(:,ia), lenTi, 1, lenia), 1, lenia, 1)./repmat(Ti.^1.5, 1, lenia, lenia);
    nuii(:, ia) = sum(nust, 3);

    % ���
    if ismember(mode, mode_i)
        varargout{1} = nuii;
    else
        varargout{2} = nuii;
    end
end

end