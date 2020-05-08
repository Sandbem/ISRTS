%%INFO: ʹ��Bessel�������Gordeyev���֡�
%%----------------------------------------------------------------------%%
% Needs: fadf.m
%%----------------------------------------------------------------------%%
% Inputs:
%   theta       - �Ա���
%   psin        - �����Գɷ���ײƵ�� [Hz]
%   alpha       - ɢ���ʸ��ų��н� [rad]
%   phi         - ����Ƶ�� [Hz]
%   N           - ���������
% Output:
%   J           - ���ֽ��
%%----------------------------------------------------------------------%%
% author: Washy [IGG]
% date: 2020/05/07
%%----------------------------------------------------------------------%%

function J = gordeyevBessel(theta,psin,alpha,phi,N)

dN = 1;
n = -N:dN:N;

eta = sin(alpha)^2/2/phi^2;
chi = besseli(n,eta,1);

n = repmat(n', 1, length(theta));
theta = repmat(theta, 2*N/dN+1, 1);

thetan = (theta - 1i*psin - n*phi)/cos(alpha);
G = sqrt(pi)/cos(alpha)*fadf(-thetan);

J = chi*G*dN;

end