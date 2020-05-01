%%INFO: ISR_mian.m���Ӻ���, ����ISR�����ס�(Mac)
%%----------------------------------------------------------------------%%
% Needs: fadf.m
%%----------------------------------------------------------------------%%
% Inputs:
%   Je          - ���ӹ�һ��Gordeyev����
%   Ji          - ���ӹ�һ��Gordeyev����
%   omge        - ���Ӷ�����Ƶ�� [rad*Hz]
%   omgi        - ���Ӷ�����Ƶ�� [rad*Hz]
%   k           - ɢ���ʸ [m^-1]
%   ne          - �������ܶ� [m^-3]
%   he          - ���ӵ°ݰ뾶 [m]
%   eta         - �����Ӵ��������ϵ��
%   mu          - eta * Tr
%   vTe         - �������ٶ� [m/s]
%   vTi         - �������ٶ� [m/s]
%   psien       - ����-���Թ�һ����ײƵ��
%   psiin       - ����-���Թ�һ����ײƵ��
% Outputs:
%   sigma       - ΢��ɢ����� [s/m]
%%----------------------------------------------------------------------%%
%�ο�����
% [1] Swartz, W. E. & Farley, D. T. (1979). Theory of Incoherent-Scattering
%   of Radio-Waves by a Plasma .5. Use of the Nyquist Theorem in General 
%   Quasi-Equilibrium Situations. Journal of Geophysical Research-Space 
%   Physics, 1979, 84, 1930-932. doi:10.1029/JA084iA05p01930
% [2] Dougherty, J. P. & Farley, D. T. (1963). A Theory of Incoherent 
%   Scattering of Radio Waves by a Plasma .3. Scattering in a Partly 
%   Ionized Gas. Journal of Geophysical Research, 1963, 68, 5473-5486.
%   doi:10.1029/jz068i019p05473
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/10
%%----------------------------------------------------------------------%%

function sigma = ISR_spectrumMac(Je,Ji,omge,omgi,k,ne,he,eta,mu,vTe,vTi,psien,psiin)

% ������Ӱ뾶re��ƽ��re2
re2 = 7.940478951426530e-30; % e^2/(4*pi*eps0*me*c^2)

% ���ӹ�һ���Ա���
thee = omge/(k*vTe); % [2] eq(7)
% ���ӹ�һ������
ye = (1i+(thee-1i*psien).*Je)./(1-psien*Je); % [2] eq(5)
% ���˶�����ĵ�����ɢ����
specTe = computeSpecT(ye,omge,psien,k,vTe,1);

% ���˶������������ɢ����
row    = length(eta);  % ����
column = length(omgi); % ����
specTi = zeros(1, column);
muyi   = complex(zeros(1, column));

for irow = 1:row
    % ���ӹ�һ���Ա���
    thei = omgi/(k*vTi(irow)); % [2] eq(7)
    % ���ӹ�һ������
    yi = (1i+(thei-1i*psiin(irow)).*Ji(irow,:))./(1-psiin(irow)*Ji(irow,:)); % [2] eq(5)
    
    specTi = specTi + computeSpecT(yi,omgi,psiin(irow),k,vTi(irow),eta(irow));
    muyi   = muyi + mu(irow)*yi;
end

den   = abs(ye+muyi+1i*k^2*he^2).^2;
spece = abs(muyi+1i*k^2*he^2).^2./den.*specTe;
speci = abs(ye).^2./den.*specTi;
% sigma = ne*re2*(spece + speci); % [1] eq(20)
sigma = 2*pi*(spece + speci); % [1] eq(20)

end

%%----------------------------------------------------------------------%%
%%INFO: ISR_spectrumMac.m���ں��Ӻ����������������˶������ɢ���ף�������
%   ���ĵ�λ��ΪNaN�����⡣
%%----------------------------------------------------------------------%%
% Needs: fadf.m
%%----------------------------------------------------------------------%%
% Inputs:
%   y           - �����ٵ���
%   omg         - ������Ƶ�� [rad*Hz]
%   psin        - �����Գɷ���ײƵ�� [Hz]
%   k           - ɢ���ʸ [m^-1]
%   vT          - ���ٶ� [m/s]
%   eta         - �����Ӵ��������ϵ��
% Outputs:
%   specT       - ���˶��������ɢ���� [s]
%%----------------------------------------------------------------------%%
% author: Washy[IGG]
% date: 2019/10/10
%%----------------------------------------------------------------------%%

function specT = computeSpecT(y,omg,psin,k,vT,eta)

specT = eta/pi*real(y)./omg; % [1] eq(20)

% ���¼���NaNλ�õ���ֵ
if ~all(omg)
    imagG = sqrt(pi)*fadf(1i*psin);
    specT(omg==0) = eta/(pi*k*vT)*imagG/(1-psin*imagG);
end

end