%%INFO: ������ACFs��
%%----------------------------------------------------------------------%%
% Input:
%   t           - �Ա���
%   funPar      - �����������
% Outputs:
%   ft          - ������
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Kudeki, E. , & Milla, M. . (2006). Incoherent scatter spectrum 
%   theory for modes propagating perpendicular to the geomagnetic field. 
%   journal of geophysical research space physics, 111(A6), -.
%   Doi: 10.1029/2005JA011546.
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/30
%%----------------------------------------------------------------------%%

function ft = SingleParticleACF(t, funPar)

psin     = funPar(1); % ��������ײƵ��
phi      = funPar(2); % ����Ƶ��
alpha    = funPar(3); % ��ų��н�
psi_par  = funPar(4); % ƽ���ڴų��Ŀ�����ײƵ��
psi_perp = funPar(5); % ��ֱ�ڴų��Ŀ�����ײƵ��

ACF_BGK = exp(-psin*t);

if psi_par+psi_perp < 1e-6 %�޿�����ײ
    ACF_par  = exp(-0.25*cos(alpha)^2*t.^2);
    ACF_perp = exp(-sin(alpha)^2*sin(0.5*phi*t).^2/phi^2);
else %�п�����ײ
    ACF_par  = exp(-0.5*cos(alpha)^2*(psi_par*t-1+exp(-psi_par*t))/psi_par^2);
    gamma    = atan(psi_perp/phi); %[1] eq(12)
    ACF_perp = exp(-0.5*sin(alpha)^2*(cos(2*gamma)+psi_perp*t-exp(-psi_perp*t).*cos(phi*t-2*gamma))/(psi_perp^2+phi^2));
end

ft = ACF_BGK.*ACF_par.*ACF_perp; %[1] eq(10)

end