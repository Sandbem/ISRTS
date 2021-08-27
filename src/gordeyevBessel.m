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

if nargin < 5
    N = 8;
end

%%
Nm = 4096;
epsJN = 1e-6;

eta = 0.5*sin(alpha)^2/phi^2;

while 1
    JN = abs(besseli(N,eta,1));
%     JN = abs(besseli(N,eta,1)*sqrt(pi)/cos(alpha)*fadf(N*phi));
    if JN < epsJN || N >= Nm
        break
    end
    N = 2*N;
end

%%
n = -N:N;
chi = besseli(n,eta,1);

n = repmat(n', 1, length(theta));
theta = repmat(theta, 2*N+1, 1);

thetan = (theta - 1i*psin - n*phi)/cos(alpha);
G = sqrt(pi)/cos(alpha)*fadf(-thetan);

J = chi*G;

end