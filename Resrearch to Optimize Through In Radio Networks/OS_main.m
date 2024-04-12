%function output = OS_main(mu_p,mu_p_1,mu_p_2,mu_p_3, h1,h2,h3) % 0<mui<0.5, 0<hi<40
clc;
clear all;

N = 4;
S = [0 0 0];
U1 = [1 2 0];
U2 = [4 2 0];
U3 = [6 2 0];
U4 = [8 2 0];

mu_p=0.1;
mu_p_1=0.2;
mu_p_2=0.3;
mu_p_3=0.4;
mu_p_4=0.5;

%h1= 0:4:40;
%h2=0:4:40;
%h3=0:4:40;
%h4=0:4:40;

h1=10;
h2=15;
h3=20;
h4=25;

P = [1 1 0];
G = [8 8 0];
Dq = [20 25 0];
Dp = [10 15 0];
% phi = [4.8860 9.6177 12.0870];
 %psi = [0.4290 0.1581 0.1139];
phi = 9.6177;
psi = 0.1581;

d_GP = 5;
d_SP = 5;
d_GDp = 8;
d_GDq = 8;

% theta_LoS = [0.1 1.0 1.6];
 %theta_NLoS = [21 20 23];
theta_LoS = 1;
theta_NLoS = 20;

m = 2;
mu_q = 1 - mu_p;
mu_q_1 = 1 - mu_p_1;
mu_q_2 = 1 - mu_p_2;
mu_q_3 = 1 - mu_p_3;
mu_q_4 = 1 -mu_p_4;    


c = 3*10^8;
fc = 2*10^9;
pi = 3.14;

% Omega_e = [5 7 9];
Omega_e = 3;
PSdB = 0:2:20;
PPdB = 0:2:20;
PU1dB = 0:2:20;
PU2dB = 0:2:20;
PU3dB = 0:2:20;
PU4dB = 0:2:20;

%PSdB = 10;
%PPdB = 10;
%PU1dB = 10;
%PU2dB = 10;
%PU3dB = 10;
%PU4dB = 10;

PS = 10.^(PSdB./10);
PP = 10.^(PPdB./10);
PU1 = 10.^(PU1dB./10);
PU2 = 10.^(PU2dB./10);
PU3 = 10.^(PU3dB./10);
PU4 = 10.^(PU4dB./10);

    d_SU1 = sqrt((U1(1))^2 + (U1(2))^2+(h1).^2);
    d_U1U2 = sqrt((U1(1)-U2(1))^2+(U1(2)-U2(2)).^2+(h1-h2).^2);
    d_U2U3 = sqrt((U2(1)-U3(1))^2+(U2(2)-U3(2)).^2+(h2-h3).^2);
    d_U3U4 = sqrt((U3(1)-U4(1))^2+(U3(2)-U4(2)).^2+(h3-h4).^2);
    d_U4Dp = sqrt((U4(1)-Dp(1))^2+(U4(2)-Dp(2)).^2+(h4).^2);
    d_U4Dq = sqrt((U4(1)-Dq(1))^2+(U4(2)-Dq(2)).^2+(h4).^2);
    d_U1P = sqrt((U1(1)-P(1))^2+(U1(2)-P(2)).^2+(h1).^2);
    d_GU4 = sqrt((G(1)-U4(1))^2+(G(2)-U4(2)).^2+(h4).^2);
   
%   d_SU1 = 3;
%   d_U1U2 = 3;
%   d_U2U3 = 3;
%   d_U3U4 = 3;
%   d_U4Dp = 3;
%   d_U4Dq = 3;
%   d_U1P = 3;
%   d_GU4 = 3;
   
nol = 10^5;

R_Oq = 2;
R_Op = 2;
R_O_P = 11;
W =10^7;


beta_SU1 = atan(h1./d_SU1);
beta_U4Dq = atan(h3./d_U4Dq);
beta_U4Dp = atan(h3./d_U4Dp);
beta_GU4 = atan(h3./d_GU4);

K_LoS = theta_LoS.*(c/(4*pi*fc))^-1;
K_NLoS = theta_NLoS.*(c/(4*pi*fc))^-1;

k_SU1 = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_SU1.*180./pi + psi.*phi));

L_SU1 = (k_SU1).*(d_SU1.^2);

k_U4Dq = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U4Dq.*180./pi + psi.*phi));

L_U4Dq = (k_U4Dq).*(d_U4Dq.^2);

k_U4Dp = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U4Dp.*180./pi + psi.*phi));

L_U4Dp = (k_U4Dp).*(d_U4Dp.^2);

L_U1U2 = 4*pi*fc/c*(d_U1U2.^2);
L_U2U3 = 4*pi*fc/c*(d_U2U3.^2);
L_U3U4 = 4*pi*fc/c*(d_U3U4.^2);

k_GU4 = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_GU4.*180./pi + psi.*phi));

L_GU4 = (k_GU4).*(d_GU4.^2);


Xi_Oq = (2^((N+1)*R_Oq/W))-1;
Xi_Op = (2^((N+1)*R_Op/W))-1;
Xi_O_P = (2^((N+1)*R_O_P/W))-1;
%Xi_O_P = 10^(-3);
% Best Relay Scheme
        
        gSU1 = max(gamrnd(m,1/m,1,1,nol),[],1);
        gU1U2= max(gamrnd(m,1/m,1,1,nol),[],1);
        gU2U3 = max(gamrnd(m,1/m,1,1,nol),[],1);
        gU3U2 = max(gamrnd(m,1/m,1,1,nol),[],1);
        gU3U4 = max(gamrnd(m,1/m,1,1,nol),[],1);
        gU4Dp = max(gamrnd(m,1/m,1,1,nol),[],1);
        gU4Dq = max(gamrnd(m,1/m,1,1,nol),[],1);
        gGU4 = max(gamrnd(m,1/m,1,1,nol),[],1);
        
        gGDp = max(gamrnd(1,1,1,1,nol),[],1);
        gGDq = max(gamrnd(1,1,1,1,nol),[],1);
        
        gGP = max(gamrnd(1,1,1,1,nol),[],1);
        gSP = max(gamrnd(1,1,1,1,nol),[],1);
        
        gU1P = max(gamrnd(m,1/m,1,1,nol),[],1);
        
%        slot 1
        gamma_U1_p = mu_p.*PS.*gSU1./(L_SU1.*(PS*Omega_e+1));
        gamma_U1_q = mu_q.*PS.*gSU1./(L_SU1.*((mu_p.*PS.*gSU1./L_SU1)+PS*Omega_e+1));
        
        gamma_1_P = PP.*gGP./((d_GP^2).*((PS.*gSP./(d_SP^2))+(PS+PP)*Omega_e)+1);
        
%        slot 2
        gamma_U2_p = mu_p_1.*PU1.*gU1U2./(L_U1U2.*(PU1*Omega_e+1));
        gamma_U2_q = mu_q_1.*PU1.*gU1U2./(L_U1U2.*((mu_p_1.*PU1.*gU1U2./L_U1U2)+PU1*Omega_e+1));
        
        gamma_2_P = PP.*gGP./((d_GP.^2).*((PU1.*gU1P./(d_U1P.^2))+(PU1+PP)*Omega_e)+1);
        
%        slot 3
        
        gamma_U3_p =  mu_p_2.*PU2.*gU2U3./(L_U2U3.*(PU2*Omega_e+1));
        gamma_U3_q = mu_q_2.*PU2.*gU2U3./(L_U2U3.*((mu_p_2.*PU2.*gU2U3./L_U1U2)+PU2*Omega_e+1));
        
        gamma_3_P = PP.*gGP./((d_GP.^2).*((PU2.*gU1P./(d_U1P.^2))+(PU1+PP)*Omega_e)+1);
        
  %      slot 4
         gamma_U4_p = mu_p_3.*PU3.*gU3U4./(L_U3U4.*((PP.*gGU4./L_GU4)+(PU3+PP)*Omega_e+1));
        gamma_U4_q = mu_q_3.*PU3.*gU3U4./(L_U3U4.*((mu_p_3.*PU3.*gU3U4./L_U3U4)+(PP.*gGU4./L_GU4)+(PU3+PP)+PU3*Omega_e+1));
        
        gamma_4_P = PP.*gGP./((d_GP^2).*(PP*Omega_e+1));
   %    Slot 5
        gamma_D_p = mu_p_4.*PU4.*gU4Dp./(L_U4Dp.*((PP.*gGDp./d_GDp^2)+(PP+PU4)*Omega_e+1));
        gamma_D_q = mu_q_4.*PU4.*gU4Dq./(L_U4Dq.*((mu_p_4.*PU4.*gU4Dq./L_U4Dq)+(PP.*gGDq./d_GDq^2)+(PP+PU4)*Omega_e+1));
        
        gamma_5_P = PP.*gGP./((d_GP^2).*(PP*Omega_e+1));
        
        
        gamma_E2E_p1 = min(gamma_U1_p,gamma_U2_p);
        gamma_E2E_p2 = min(gamma_U3_p,gamma_U4_p);
        gamma_E2E_p3 = min(gamma_D_p,gamma_E2E_p2);
        gamma_E2E_p = min(gamma_E2E_p1,gamma_E2E_p3);
        
        gamma_E2E_q1 = min(gamma_U1_q,gamma_U2_q);
        gamma_E2E_q2 = min(gamma_U3_q,gamma_U4_q);
        gamma_E2E_q3 = min(gamma_E2E_q2,gamma_D_q);
        gamma_E2E_q = min(gamma_E2E_q1,gamma_E2E_q3);

        gamma_P1 = min(gamma_1_P,gamma_2_P);
        gamma_P2 = min(gamma_3_P,gamma_4_P);
        gamma_P3 = min(gamma_P2,gamma_5_P);
        gamma_P = min(gamma_P1,gamma_P3);
    countD = 0;
    for  j = 1:nol
        if gamma_E2E_p(j) <= Xi_Op || gamma_E2E_q(j)<=Xi_Oq
            countD = countD +1;
            
        end
    end
   
    O_D = countD/nol;
%    throughput of Secondary network
    Th_D = (1-O_D)*R_Oq; 
 %   minimun Th_D
      %mu_p,mu_p_1,mu_p_2, mu_p_3
         countP = 0; 
         O_P=0;
       for j = 1:length(PPdB)
            for  i = 1:nol
                if gamma_P(i) <= Xi_O_P 
                    countP = countP +1;
                    
                    O_P(j) = countP/nol;
                   
                    %countP
                   % O_P
                end
            end
                     % countP
%           outage probability of Primary network
        
           
  %   Constraint O_P < epsilon --> PS, PU1,
    end
% output = Th_q;
for i = 1:length(PSdB)
  %O_PO1(i) =  tinhtoan(gamma_P,PP(i),Omega_e,W,d_GP);
    O_PO1_MP(i) = mophong(gamma_P,m,U1,h1,N,S,P,G,PP(i),Omega_e,R_O_P,W,nol);
   % O_PO1_MP
end
%subplot(2,1,1);
%plot(O_PO1_MP,'marker','o','markersize',5,'markerfacecolor','r');
%subplot(2,1,2);
%plot(O_P,'marker','s','markersize',5,'markerfacecolor','r');


%semilogy (PPdB,O_P,'r');
%legend('O_PO1_MP', 'O_P');
%ylabel('Xac suat dung');
%xlabel('Cong suat phat');

semilogy (PSdB,O_P,'r',PSdB,O_P,'ro');
ylabel('Xac suat dung');
xlabel('Cong suat');
grid on;
%grid on
% output = O_P;
%output = Th_D;
%output = d_GU4;

    