function OP_CGA()

% clc;
clear all;
close all
%^=====================================================
count=0;
r=20;%number of generation
psize=20;
% N=[10,50];% population size
n=31; % number of optimization variables 
load data1 data1;
block=1;
x=data1(1:psize,:)';
% x=u';
bound    =[     1  1   1 0.01 0.1 0.1 0.1 0.1 0.1 0.1 1 2 0 4 2 0 6 2 0 8  3 0 10 4 0 0  0  0  0  0  3  ;
            100 100 100 0.99 0.5 0.5 0.5 0.5 0.5 0.5 4 5 0 7 6 0 9 6 0 11 7 0 13 8 0 40 40 40 40 40 5 ];
bound = bound';
mutrate=0.03; %  mutation rate 
sigma=0.02; % mutation deviation
selection=0.0;  %fraction of population kept 
crosrate=0.7;%crossover rate
% rng(10)% random seed
% for jj=1:length(N)
%     x1 = rand(1,N(jj))*3;
%     x2 = rand(1,N(jj))*3;
%     p=[x1;x2];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PTOP =[];
    for j=1:r %generation num
        for i=1:psize
            x(:,i) = max(x(:,i),bound(:,1));
            x(:,i) = min(x(:,i),bound(:,2));
            x(end,i) = round(x(end,i));
            OP_thea(i) = Ob_x(x(:,i));
            y=OP_thea;
        end
        
        mean_mat3(j) = mean(OP_thea);
        
        OP_thea = OP_thea - min(OP_thea);
        [idd,order] = sort(OP_thea);
        best_cand(:,j) = x(:,order(psize));
        pop_select= select(x,OP_thea); %select population by roullet   
        
  %%%% Crossover      
        pp = pop_select;
        tt = 1;
        while tt<psize*crosrate
            beta = rand(1);
            pp(:,tt)=beta* pop_select(:,tt)   + (1-beta)*pop_select(:,tt+1);
            pp(:,tt+1)=beta* pop_select(:,tt+1)   + (1-beta)*pop_select(:,tt);
            tt = tt+2;
        end
        pop_select = pp;
 %%% Mutation
        for kk = 1:ceil(mutrate*psize)
            picked_mute = ceil(rand(1)*psize);
            pop_select(:,picked_mute) = pop_select(:,picked_mute) +  sigma*randn(n,1);
        end
        x = pop_select;
    
    end
    figure(1)
    plot(1:r, mean_mat3,'ro-');
    hold on;
    best_cand = best_cand;
    save mean_mat3 mean_mat3;   
    best_of_best = best_cand(:,end);
    disp('==================================');
    disp(Ob_x(best_of_best))


function y = Ob_x (x)
PP = x(1);
PS = x(2);
PB = x(3);
tau = x(4);
mu_p = x(5);
mu_p_1 = x(6);
mu_p_2 = x(7);
mu_p_3 = x(8);
mu_p_4 = x(9);
% = x(10);
U1 = x(11:13);
% U1(1)=x(11);
% U1(2)=x(12);
% U1(3)=x(13);
U2 = x(14:16);
% U2(1)=x(14);
% U2(2)=x(15);
% U2(3)=x(16);
U3 = x(17:19);
% U3(1)= x(17);
% U3(2)= x(18);
% U3(3)= x(19);
U4 = x(20:22);
% U4(1) = x(20);
% U4(2) = x(21);
% U4(3) = x(22);
%U5 = x(23:25);
% U5(1) = x(23);
% U5(2) = x(24);
% U5(3) = x(25);
h1 = x(26);
h2 = x(27);
h3 = x(28);
h4 =x(29);
% = x(30);
N = x(31);
y = Objective_F(PP,PS,PB,tau,mu_p,mu_p_1,mu_p_2,mu_p_3,mu_p_4,,U1, U2, U3, U4, U5, h1,h2,h3,h4,N);
%    
%    
% x = [PP,PS,PB,tau,mu_p,mu_p_1,mu_p_2,mu_p_3,mu_p_4,,U1, U2, U3, U4, U5, h1,h2,h3,h4,,N];
% y = Objective_F(PP,PS,PB,tau,mu_p,mu_p_1,mu_p_2,mu_p_3,mu_p_4,,...
%     U1(1),U1(2),U1(3),U2(1),U2(1),U2(1), U3(1),U3(2),U3(3), U4(1),U4(2),U4(3), U5(1),U5(2),U5(3),h1,h2,h3,h4,,N);
% function [y_out] =  Objective_F(PP,PS,PB,tau,mu_p,mu_p_1,mu_p_2,mu_p_3,mu_p_4,,U1(1),U1(2),U1(3), U2(1),U2(1),U2(1),U3(1),U3(2),U3(3), U4(1),U4(2),U4(3), U5(1),U5(2),U5(3),h1,h2,h3,h4,,N)
% % function y_out = Objective_F(PP,PS,PB,tau,mu_p,mu_p_1,mu_p_2,mu_p_3,mu_p_4,,...
% % %        U1, U2, U3, U4, U5, h1,h2,h3,h4,,N)

function y_out = Objective_F(PP,PS,PB,tau,mu_p,mu_p_1,mu_p_2,mu_p_3,mu_p_4,,U1,U2,U3,U4,U5,h1,h2,h3,h4,,N)
S = [0 0 0];
P = [1 1 0];
G = [8 8 0];
Dq = [15 25 0];
Dp = [10 15 0];
E = [25 15 0];

    d_SU1 = sqrt((U1(1))^2 + (U1(2))^2+(h1).^2);
    d_U1U2 = sqrt((U1(1)-U2(1))^2+(U1(2)-U2(2)).^2+(h1-h2).^2);
    d_U2U3 = sqrt((U2(1)-U3(1))^2+(U2(2)-U3(2)).^2+(h2-h3).^2);
    d_U3U4 = sqrt((U3(1)-U4(1))^2+(U3(2)-U4(2)).^2+(h3-h4).^2);
   %  = sqrt((U4(1)-U5(1))^2+(U4(2)-U5(2)).^2+(h4-).^2);
    
    d_U3Dp = sqrt((U3(1)-Dp(1))^2+(U3(2)-Dp(2)).^2+(h3).^2);
    
    d_U3Dq = sqrt((U3(1)-Dq(1))^2+(U3(2)-Dq(2)).^2+(h3).^2);
    
    d_U4Dp = sqrt((U4(1)-Dp(1))^2+(U4(2)-Dp(2)).^2+(h4).^2);
    d_U4Dq = sqrt((U4(1)-Dq(1))^2+(U4(2)-Dq(2)).^2+(h4).^2);
    
   %  = sqrt((U5(1)-Dp(1))^2+(U5(2)-Dp(2)).^2+().^2);
   %  = sqrt((U5(1)-Dq(1))^2+(U5(2)-Dq(2)).^2+().^2);
    
    d_U1P = sqrt((U1(1)-P(1))^2+(U1(2)-P(2)).^2+(h1).^2);
    d_GU3 = sqrt((G(1)-U3(1))^2+(G(2)-U3(2)).^2+(h3).^2);

    d_GU4 = sqrt((G(1)-U4(1))^2+(G(2)-U4(2)).^2+(h4).^2);
   %  = sqrt((G(1)-U5(1))^2+(G(2)-U5(2)).^2+().^2);
    
d_GP = sqrt((G(1)-P(1))^2+(G(2)-P(2))^2);
d_SP = sqrt((S(1)-P(1))^2+(S(2)-P(2))^2);
d_GDp = sqrt((G(1)-Dp(1))^2+(G(2)-Dp(2))^2);
d_GDq = sqrt((G(1)-Dq(1))^2+(G(2)-Dq(2))^2);

%%%%%Distance EAV%%%%%%%%%%%%%%%%
d_U3E = sqrt((U3(1)-E(1))^2+(U3(2)-E(2)).^2+(h3).^2);
d_U4E = sqrt((U4(1)-E(1))^2+(U4(2)-E(2)).^2+(h4).^2);
% = sqrt((U5(1)-E(1))^2+(U5(2)-E(2)).^2+().^2);
d_GE = sqrt((G(1)-E(1))^2+(G(2)-E(2))^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
phi = [4.8860 9.6177 12.0870]; phi = phi(2);
psi = [0.4290 0.1581 0.1139]; psi = psi(2);
theta_LoS = [0.1 1.0 1.6]; theta_LoS = theta_LoS(2);
theta_NLoS = [21 20 23];theta_NLoS = theta_NLoS(2);
m = 2;


c = 3*10^8;
fc = 2*10^9;
pi = 3.14;

Omega_e = [3 5 7]; Omega_e = Omega_e(2);
Omega_PB = 10;
mu_q = 1 - mu_p;
mu_q_1 = 1 - mu_p_1;
mu_q_2 = 1 - mu_p_2;
mu_q_3 = 1 - mu_p_3;
mu_q_4 = 1 - mu_p_4;
 = 1 - ;
nol = 10^5;
delta = 0.85;

R_Oq = 10;
R_Op = 10;
R_O_P = 1000;
R_E = 10;
W =10^8;

V=[2 3 4] %%%%%-------%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
        y_out =OP_the(PP,PS,PB,tau,delta,mu_p,mu_p_1,mu_p_2,mu_p_3,mu_p_4,,...
        mu_q,mu_q_1,mu_q_2,mu_q_3,mu_q_4,,d_SU1,d_U1U2,d_U2U3,d_U3Dp,d_U3Dq,d_U1P,d_GU3,...
        d_U3U4,,d_U4Dp,d_U4Dq,,,d_GU4,,d_U3E,d_U4E,,d_GE,...
        d_GP,d_SP,d_GDp,d_GDq,W,fc,c,pi,phi,psi,theta_LoS,...
        theta_NLoS,R_Oq,R_Op,R_O_P,R_E,h1,h2,h3,h4,,N,m,Omega_e,Omega_PB,nol);
   
% end
function output = OP_the(PP,PS,PB,tau,delta,mu_p,mu_p_1,mu_p_2,mu_p_3,mu_p_4,,...
        mu_q,mu_q_1,mu_q_2,mu_q_3,mu_q_4,,d_SU1,d_U1U2,d_U2U3,d_U3Dp,d_U3Dq,d_U1P,d_GU3,...
        d_U3U4,,d_U4Dp,d_U4Dq,,,d_GU4,,d_U3E,d_U4E,,d_GE,...
        d_GP,d_SP,d_GDp,d_GDq,W,fc,c,pi,phi,psi,theta_LoS,...
        theta_NLoS,R_Oq,R_Op,R_O_P,R_E,h1,h2,h3,h4,,N,m,Omega_e,Omega_PB,nol)

V = 2; % ADD MORE ANTTENNA
I_P = 5; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_SU1 = atan(h1./d_SU1);
beta_U1P = atan(h1./d_U1P);

beta_U3Dq = atan(h3./d_U3Dq);
beta_U3Dp = atan(h3./d_U3Dp);
beta_GU3 = atan(h3./d_GU3);

beta_U4Dq = atan(h4./d_U4Dq);
beta_U4Dp = atan(h4./d_U4Dp);
beta_GU4 = atan(h4./d_GU4);

%beta_U5Dq = atan(./);
%beta_U5Dp = atan(./);
%beta_GU5 = atan(./);

%%%%%%%%%%%%%%%%EAV%%%%%%%%%%%%%
beta_U3E = atan(h3./d_U3E);
beta_U4E = atan(h4./d_U4E);
%beta_U5E = atan(./d_U4E);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_LoS = theta_LoS.*(c/(4*pi*fc))^-1;
K_NLoS = theta_NLoS.*(c/(4*pi*fc))^-1;

k_SU1 = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_SU1.*180./pi + psi.*phi));

L_SU1 = (k_SU1)*(d_SU1^2);

k_U1P = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U1P.*180./pi + psi.*phi));

L_U1P = (k_U1P)*(d_U1P^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%

k_U3Dq = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U3Dq.*180./pi + psi.*phi));

L_U3Dq = (k_U3Dq)*(d_U3Dq^2);

k_U3Dp = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U3Dp.*180./pi + psi.*phi));

L_U3Dp = (k_U3Dp)*(d_U3Dp^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%
k_U4Dq = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U4Dq.*180./pi + psi.*phi));

L_U4Dq = (k_U4Dq)*(d_U4Dq^2);

k_U4Dp = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U4Dp.*180./pi + psi.*phi));

L_U4Dp = (k_U4Dp)*(d_U4Dp^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%k_U5Dq = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U5Dq.*180./pi + psi.*phi));

%L_U5Dq = (k_U5Dq)*(^2);

%k_U5Dp = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U5Dp.*180./pi + psi.*phi));

%L_U5Dp = (k_U5Dp)*(^2);
%%%%%%%%%%%%%%%%%%%%

L_U1U2 = 4*pi*fc/c*(d_U1U2^2);
L_U2U3 = 4*pi*fc/c*(d_U2U3^2);
L_U3U4 = 4*pi*fc/c*(d_U3U4^2);
%L_U4U5 = 4*pi*fc/c*(^2);

k_GU3 = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_GU3.*180./pi + psi.*phi));
k_GU4 = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_GU4.*180./pi + psi.*phi));
%k_GU5 = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_GU5.*180./pi + psi.*phi));


L_GU3 = (k_GU3)*(d_GU3^2);
L_GU4 = (k_GU4)*(d_GU4^2);
%L_GU5 = (k_GU5)*(^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_U3E = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U3E.*180./pi + psi.*phi));
L_U3E = (k_U3E)*(d_U3E^2);
k_U4E = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U4E.*180./pi + psi.*phi));
L_U4E = (k_U4E)*(d_U4E^2);
%k_U5E = K_NLoS + (K_LoS-K_NLoS)./(1+phi.*exp(-psi.*beta_U5E.*180./pi + psi.*phi));
%L_U5E = (k_U5E)*(^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xi_Oq = (2^((N+1)*R_Oq/((1-tau)*W)))-1;
Xi_Op = (2^((N+1)*R_Op/((1-tau)*W)))-1;
Xi_O_P = (2^((N+1)*R_O_P/((1-tau)*W)))-1;
Xi_E = (2^((N+1)*R_E/((1-tau)*W)))-1;

Omega_GP =1;
Omega_SP =1;
Omega_U1P =1;
Omega_SU1 = 1;
Omega_U1U2 = 1;
Omega_GU3 = 1;
Omega_GU4 = 1;
%Omega_GU5 = 1;

Omega_U2U3 = 1;
Omega_U3U4 = 1;
%Omega_U4U5 = 1;

Omega_GDp = 1;
Omega_U3Dp = 1;
Omega_U4Dp = 1;
%Omega_U5Dp = 1;

Omega_GDq= 1;
Omega_U3Dq = 1;
Omega_U4Dq = 1;
%Omega_U5Dq = 1;
Omega_U3E =1;
Omega_U4E =1;
Omega_U5E =1;
Omega_GE =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if N==3
    %%%%%%%%%%EH at UAV%%%%%%%%%%%
    PU1 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    PU2 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    PU3 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    
    
    %%%%%%%%%%%% OP phase 1 %%%%%%%%%%%%%%%%%%
    tu_1_P = exp(-Xi_O_P*d_GP^2*((PS+PP)*Omega_e+1)/(PP*Omega_GP));
    mau_1_P = Omega_SP*((Xi_O_P*d_GP^2*PS/(d_SP^2*PP*Omega_GP))+1/Omega_SP);
    O_1_P = 1-tu_1_P/mau_1_P;
    
    %%%%%%%%%%%% OP phase 2 %%%%%%%%%%%%%%%%%%
    tu_2_P = ((m/Omega_U1P)^m)*exp(-Xi_O_P*d_GP^2*((PU1+PP)*Omega_e+1)/(PP*Omega_GP))*factorial(m-1);
    mau_2_P = gamma(m)*((Xi_O_P*d_GP^2*PU1/(L_U1P*PP*Omega_GP))+m/Omega_U1P)^2;
    O_2_P = 1-tu_2_P/mau_2_P;
    
    %%%%%%%%%%%% OP phase n %%%%%%%%%%%%%%%%%%
    tu_3_P = exp(-Xi_O_P*d_GP^2*(PP*Omega_e+1)/(PP*Omega_GP));
    O_3_P = 1-tu_3_P;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Throughput of the secondary network Dp
    
    P1_p = 0;
    P2_p = 0;
    P4_p_temp = 0;
    P5_p_temp =0;
    P4_p_tempE = 0;
    %P5_p_tempE =0;
    for j = 0 : m-1
        A1 = m*Xi_Op*L_SU1*(PS*Omega_e + I_P +1)/(Omega_SU1*mu_p*PS);
        P1_p =P1_p + (A1^j)*exp(-A1)/factorial(j);
        A2 =  m*Xi_Op*L_U1U2*(PU1*Omega_e + I_P +1)/(Omega_U1U2*mu_p_1*PU1);
        P2_p =P2_p + (A2^j)*exp(-A2)/factorial(j);
        %%%%%%%%%%%%%%%
        Temp_3_p = 0;
        for k = 0:j
            A3_tu = nchoosek(j,k)*(((PP+PU2)*Omega_e + I_P +1)^(j-k))*((PP/L_GU3)^k)*factorial(k+m-1);
            A3_mau = gamma(m)*(((m*Xi_Op*L_U2U3*PP/(Omega_U2U3*mu_p_2*PU2*L_GU3))+m/Omega_GU3)^(k+m));
            Temp_3_p = Temp_3_p + A3_tu/A3_mau;
        end
        P4_p_temp = P4_p_temp + (Temp_3_p*(m*Xi_Op*L_U2U3/(Omega_U2U3*mu_p_2*PU2))^j)/factorial(j);
        %%%%%%%%%%%%
       % Temp_5_p = 0;
       % for k = 0:j
         %   A5_tu = nchoosek(j,k)*(((PP+PU3)*Omega_e + I_P +1)^(j-k))*((PP/d_GDp^2)^k)*factorial(k);
         %   A5_mau = (m*Xi_Op*L_U3Dp*PP/(Omega_U3Dp*mu_p_3*PU3*d_GDp^2)+1/Omega_GDp)^(k+1);
         %   Temp_5_p = Temp_5_p+A5_tu/A5_mau;
        %end
        P5_p_temp = P5_p_temp + (Temp_5_p*(m*Xi_Op*L_U3Dp/(Omega_U3Dp*mu_p_3*PU3))^j)/factorial(j);
        %%%%%%%%%%%%%%%%%EAV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Temp_3_pE = 0;
        for k = 0:j
            A3_tu = nchoosek(j,k)*(((PP+PU2)*Omega_e + I_P +1)^(j-k))*((PP/L_GU3)^k)*factorial(k+m-1);
            A3_mau = gamma(m)*(((m*Xi_E*L_U2U3*PP/(Omega_U2U3*mu_p_2*PU2*L_GU3))+m/Omega_GU3)^(k+m));
            Temp_3_pE = Temp_3_pE + A3_tu/A3_mau;
        end
        P4_p_tempE = P4_p_tempE + (Temp_3_pE*(m*Xi_E*L_U2U3/(Omega_U2U3*mu_p_2*PU2))^j)/factorial(j);
        %%%%%%%%%%%%
        %Temp_5_pE = 0;
        %for k = 0:j
          %  A5_tu = nchoosek(j,k)*(((PP+PU3)*Omega_e + I_P +1)^(j-k))*((PP/d_GE^2)^k)*factorial(k);
           % A5_mau = (m*Xi_E*L_U3E*PP/(Omega_U3E*mu_p_3*PU3*d_GE^2)+1/Omega_GE)^(k+1);
           % Temp_5_pE = Temp_5_pE+A5_tu/A5_mau;
       % end
       % P5_p_tempE = P5_p_tempE + (Temp_5_p*(m*Xi_E*L_U3E/(Omega_U3E*mu_p_3*PU3))^j)/factorial(j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    P3_p = P4_p_temp*((m/Omega_GU3)^m)*exp(-m*Xi_Op*L_U2U3*((PP+PU2)*Omega_e + I_P +1)/(Omega_U2U3*mu_p_2*PU2));
    P4_p = P5_p_temp*exp(-(m*Xi_Op*L_U3Dp*((PP+PU3)*Omega_e + I_P +1))/(Omega_U3Dp*mu_p_3*PU3))/Omega_GDp;
    
    P3_pE = P4_p_tempE*((m/Omega_GU3)^m)*exp(-m*Xi_Op*L_U2U3*((PP+PU2)*Omega_e + I_P +1)/(Omega_U2U3*mu_p_2*PU2));
    P4_pE = P5_p_tempE*exp(-(m*Xi_Op*L_U3E*((PP+PU3)*Omega_e + I_P +1))/(Omega_U3E*mu_p_3*PU3))/Omega_GE;
    
    
    O_Dp = 1-P1_p*P2_p*P3_p*P4_p;
    Th_Dp = (1-O_Dp)*R_Op;
    
    %%%%%%%%%Leakage EAV%%%%%%%%%
    L_Dp = P1_p*P2_p*P3_pE*P4_pE;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Throughput of the secondary network Dq
    P1_q = 0;
    P2_q = 0;
    P3_q_temp = 0;
    P4_q_temp =0;
    P4_q_tempE =0;
    for j = 0 : m-1
        B1 = m*Xi_Oq*L_SU1*(PS*Omega_e+1)/(Omega_SU1*(mu_q-mu_p*Xi_Oq)*PS);
        P1_q =P1_q + (B1^j)*exp(-B1)/factorial(j);
        B2 =  m*Xi_Oq*L_U1U2*(PU1*Omega_e+1)/(Omega_U1U2*(mu_q_1-mu_p_1*Xi_Oq)*PU1);
        P2_q =P2_q + (B2^j)*exp(-B2)/factorial(j);
        %%%%%%%%%%%%%%%
        Temp_3_q = 0;
        for k = 0:j
            B3_tu = nchoosek(j,k)*(((PP+PU2)*Omega_e+1)^(j-k))*((PP/L_GU3)^k)*factorial(k+m-1);
            B3_mau = gamma(m)*(((m*Xi_Oq*L_U2U3*PP/(Omega_U2U3*(mu_q_2-mu_p_2*Xi_Oq)*PU2*L_GU3))+m/Omega_GU3)^(k+m));
            Temp_3_q = Temp_3_q + B3_tu/B3_mau;
        end
        P3_q_temp = P3_q_temp + (Temp_3_q*(m*Xi_Oq*L_U2U3/(Omega_U2U3*(mu_q_2-mu_p_2*Xi_Oq)*PU2))^j)/factorial(j);
        %%%%%%%%%%%%
        Temp_4_q = 0;
        for k = 0:j
            B4_tu = nchoosek(j,k)*(((PP+PU3)*Omega_e+1)^(j-k))*((PP/d_GDq^2)^k)*factorial(k);
            B4_mau = (m*Xi_Oq*L_U3Dq*PP/(Omega_U3Dq*(mu_q_3-mu_p_3*Xi_Oq)*PU3*d_GDq^2)+1/Omega_GDq)^(k+1);
            Temp_4_q = Temp_4_q+B4_tu/B4_mau;
        end
        P4_q_temp = P4_q_temp + (Temp_4_q*(m*Xi_Oq*L_U3Dq/(Omega_U3Dq*(mu_q_3-mu_p_3*Xi_Oq)*PU3))^j)/factorial(j);
        %%%%%%%%%%%%EAV%%%%%%%%%%%%%%
        Temp_4_qE = 0;
        for k = 0:j
            B4_tu = nchoosek(j,k)*(((PP+PU3)*Omega_e+1)^(j-k))*((PP/d_GE^2)^k)*factorial(k);
            B4_mau = (m*Xi_Oq*L_U3E*PP/(Omega_U3E*(mu_q_3-mu_p_3*Xi_E)*PU3*d_GE^2)+1/Omega_GE)^(k+1);
            Temp_4_qE = Temp_4_qE+B4_tu/B4_mau;
        end
        P4_q_tempE = P4_q_tempE + (Temp_4_qE*(m*Xi_E*L_U3E/(Omega_U3E*(mu_q_3-mu_p_3*Xi_E)*PU3))^j)/factorial(j);
        %%%%%%%%%%%%%%%%%%%%%
    end
    P3_q = P3_q_temp*((m/Omega_GU3)^m)*exp(-m*Xi_Oq*L_U2U3*((PP+PU2)*Omega_e+1)/(Omega_U2U3*(mu_q_1-mu_p_1*Xi_Oq)*PU2));
    P4_q = P4_q_temp*exp(-(m*Xi_Oq*L_U3Dq*((PP+PU3)*Omega_e+1))/(Omega_U3Dq*(mu_q_3-mu_p_3*Xi_Oq)*PU3))/Omega_GDq;
    
    P4_qE = P4_q_tempE*exp(-(m*Xi_E*L_U3E*((PP+PU3)*Omega_e+1))/(Omega_U3E*(mu_q_3-mu_p_3*Xi_E)*PU3))/Omega_GE;
    
    O_Dq = 1-P1_q*P2_q*P3_q*P4_q;
    Th_Dq = (1-O_Dq)*R_Oq;
    
    L_Dq = P1_q*P2_q*P3_q*P4_qE;
elseif N==4
    %%%%%%%%%%EH at UAV%%%%%%%%%%%
    PU1 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    PU2 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    PU3 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    PU4 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    %%%%%%%%%%%% OP phase 1 %%%%%%%%%%%%%%%%%%
    tu_1_P = exp(-Xi_O_P*d_GP^2*((PS+PP)*Omega_e+1)/(PP*Omega_GP));
    mau_1_P = Omega_SP*((Xi_O_P*d_GP^2*PS/(d_SP^2*PP*Omega_GP))+1/Omega_SP);
    O_1_P = 1-tu_1_P/mau_1_P;
    O_1_P = O_1_P^V
    
    %%%%%%%%%%%% OP phase 2 %%%%%%%%%%%%%%%%%%
    tu_2_P = ((m/Omega_U1P)^m)*exp(-Xi_O_P*d_GP^2*((PU1+PP)*Omega_e+1)/(PP*Omega_GP))*factorial(m-1);
    mau_2_P = gamma(m)*((Xi_O_P*d_GP^2*PU1/(L_U1P*PP*Omega_GP))+m/Omega_U1P)^2;
    O_2_P = 1-tu_2_P/mau_2_P;
    O_2_P = O_2_P^V
    
    %%%%%%%%%%%% OP phase n %%%%%%%%%%%%%%%%%%
    tu_3_P = exp(-Xi_O_P*d_GP^2*(PP*Omega_e+1)/(PP*Omega_GP));
    O_3_P = 1-tu_3_P;
    O_3_P = O_3_P^V
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P1_p = 0;
    P2_p = 0;
    P3_p = 0;
    P4_p_temp = 0;
    %P5_p_temp =0;
    %P5_p_tempE =0;
    for j = 0 : m-1
        A1 = m*Xi_Op*L_SU1*(PS*Omega_e + I_P +1)/(Omega_SU1*mu_p*PS);
        P1_p =P1_p + (A1^j)*exp(-A1)/factorial(j);
        
        A2 =  m*Xi_Op*L_U1U2*(PU1*Omega_e + I_P +1)/(Omega_U1U2*mu_p_1*PU1);
        P2_p =P2_p + (A2^j)*exp(-A2)/factorial(j);
        
        A3 =  m*Xi_Op*L_U2U3*(PU2*Omega_e + I_P +1)/(Omega_U2U3*mu_p_2*PU1);
        P3_p =P3_p + (A3^j)*exp(-A3)/factorial(j);
        %%%%%%%%%%%%%%%
        Temp_3_p = 0;
        for k = 0:j
            A3_tu = nchoosek(j,k)*(((PP+PU3)*Omega_e + I_P +1)^(j-k))*((PP/L_GU4)^k)*factorial(k+m-1);
            A3_mau = gamma(m)*(((m*Xi_Op*L_U3U4*PP/(Omega_U3U4*mu_p_3*PU3*L_GU4))+m/Omega_GU4)^(k+m));
            Temp_3_p = Temp_3_p + A3_tu/A3_mau;
        end
        P4_p_temp = P4_p_temp + (Temp_3_p*(m*Xi_Op*L_U3U4/(Omega_U3U4*mu_p_3*PU3))^j)/factorial(j);
        %%%%%%%%%%%%
       % Temp_5_p = 0;
        %for k = 0:j
          %  A5_tu = nchoosek(j,k)*(((PP+PU4)*Omega_e + I_P +1)^(j-k))*((PP/d_GDp^2)^k)*factorial(k);
           % A5_mau = (m*Xi_Op*L_U4Dp*PP/(Omega_U4Dp*mu_p_4*PU3*d_GDp^2)+1/Omega_GDp)^(k+1);
           % Temp_5_p = Temp_5_p+A5_tu/A5_mau;
        %end
        %P5_p_temp = P5_p_temp + (Temp_5_p*(m*Xi_Op*L_U4Dp/(Omega_U4Dp*mu_p_4*PU4))^j)/factorial(j);
        %%%%%%%%%%%%%%%%EAV%%%%%%%%%%%%%%%%%%%%%%%
       % Temp_5_pE = 0;
       % for k = 0:j
           % A5_tu = nchoosek(j,k)*(((PP+PU4)*Omega_e + I_P +1)^(j-k))*((PP/d_GE^2)^k)*factorial(k);
           % A5_mau = (m*Xi_E*L_U4E*PP/(Omega_U4E*mu_p_4*PU3*d_GE^2)+1/Omega_GE)^(k+1);
           % Temp_5_pE = Temp_5_pE+A5_tu/A5_mau;
        %end
       % P5_p_tempE = P5_p_tempE + (Temp_5_pE*(m*Xi_E*L_U4E/(Omega_U4E*mu_p_4*PU4))^j)/factorial(j);
        %%%%%%%%%%%%%%%%%%%
        
    end
    P4_p = P4_p_temp*((m/Omega_GU4)^m)*exp(-m*Xi_Op*L_U3U4*((PP+PU3)*Omega_e + I_P +1)/(Omega_U3U4*mu_p_3*PU3));
   % P5_p = P5_p_temp*exp(-(m*Xi_Op*L_U4Dp*((PP+PU4)*Omega_e + I_P +1))/(Omega_U4Dp*mu_p_4*PU4))/Omega_GDp;
    
   % P5_pE = P5_p_tempE*exp(-(m*Xi_E*L_U4E*((PP+PU4)*Omega_e + I_P +1))/(Omega_U4E*mu_p_4*PU4))/Omega_GE;
    
    O_Dp = 1-P1_p*P2_p*P3_p*P4_p*P5_p;
    Th_Dp = (1-O_Dp)*R_Op;
    
    L_Dp = P1_p*P2_p*P3_p*P4_p*P5_pE;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Throughput of the secondary network Dq
    P1_q = 0;
    P2_q = 0;
    P3_q = 0;
    P3_q_temp = 0;
    P4_q_temp =0;
    P4_q_tempE =0;
    for j = 0 : m-1
        B1 = m*Xi_Oq*L_SU1*(PS*Omega_e+1)/(Omega_SU1*(mu_q-mu_p*Xi_Oq)*PS);
        P1_q =P1_q + (B1^j)*exp(-B1)/factorial(j);
        
        B2 =  m*Xi_Oq*L_U1U2*(PU1*Omega_e+1)/(Omega_U1U2*(mu_q_1-mu_p_1*Xi_Oq)*PU1);
        P2_q =P2_q + (B2^j)*exp(-B2)/factorial(j);
        
        B3 =  m*Xi_Oq*L_U2U3*(PU2*Omega_e+1)/(Omega_U2U3*(mu_q_2-mu_p_2*Xi_Oq)*PU2);
        P3_q =P3_q + (B3^j)*exp(-B3)/factorial(j);
        %%%%%%%%%%%%%%%
        Temp_3_q = 0;
        for k = 0:j
            B3_tu = nchoosek(j,k)*(((PP+PU3)*Omega_e+1)^(j-k))*((PP/L_GU4)^k)*factorial(k+m-1);
            B3_mau = gamma(m)*(((m*Xi_Oq*L_U3U4*PP/(Omega_U2U3*(mu_q_3-mu_p_3*Xi_Oq)*PU3*L_GU4))+m/Omega_GU4)^(k+m));
            Temp_3_q = Temp_3_q + B3_tu/B3_mau;
        end
        P3_q_temp = P3_q_temp + (Temp_3_q*(m*Xi_Oq*L_U3U4/(Omega_U3U4*(mu_q_3-mu_p_3*Xi_Oq)*PU3))^j)/factorial(j);
        %%%%%%%%%%%%
        Temp_4_q = 0;
        for k = 0:j
            B4_tu = nchoosek(j,k)*(((PP+PU4)*Omega_e+1)^(j-k))*((PP/d_GDq^2)^k)*factorial(k);
            B4_mau = (m*Xi_Oq*L_U4Dq*PP/(Omega_U4Dq*(mu_q_4-mu_p_4*Xi_Oq)*PU4*d_GDq^2)+1/Omega_GDq)^(k+1);
            Temp_4_q = Temp_4_q+B4_tu/B4_mau;
        end
        P4_q_temp = P4_q_temp + (Temp_4_q*(m*Xi_Oq*L_U4Dq/(Omega_U4Dq*(mu_q_4-mu_p_4*Xi_Oq)*PU4))^j)/factorial(j);
        
        %%%%%%%%%%%%EAV%%%%%%%%%%%%
        Temp_4_qE = 0;
        for k = 0:j
            B4_tu = nchoosek(j,k)*(((PP+PU4)*Omega_e+1)^(j-k))*((PP/d_GE^2)^k)*factorial(k);
            B4_mau = (m*Xi_E*L_U4E*PP/(Omega_U4E*(mu_q_4-mu_p_4*Xi_E)*PU4*d_GE^2)+1/Omega_GE)^(k+1);
            Temp_4_qE = Temp_4_qE+B4_tu/B4_mau;
        end
        P4_q_tempE = P4_q_tempE + (Temp_4_qE*(m*Xi_E*L_U4E/(Omega_U4E*(mu_q_4-mu_p_4*Xi_E)*PU4))^j)/factorial(j);
        %%%%%%%%%%%%%%%%%%%%%
        
    end
    P4_q = P3_q_temp*((m/Omega_GU4)^m)*exp(-m*Xi_Oq*L_U3U4*((PP+PU3)*Omega_e+1)/(Omega_U3U4*(mu_q_4-mu_p_4*Xi_Oq)*PU3));
  %  P5_q = P4_q_temp*exp(-(m*Xi_Oq*L_U4Dq*((PP+PU4)*Omega_e+1))/(Omega_U4Dq*(mu_q_4-mu_p_4*Xi_Oq)*PU4))/Omega_GDq;
    
  %  P5_qE = P4_q_tempE*exp(-(m*Xi_E*L_U4E*((PP+PU4)*Omega_e+1))/(Omega_U4E*(mu_q_4-mu_p_4*Xi_E)*PU4))/Omega_GE;
    
    O_Dq = 1-P1_q*P2_q*P3_q*P4_q*P5_q;
    Th_Dq = (1-O_Dq)*R_Oq;
    
    L_Dq = P1_q*P2_q*P3_q*P4_q*P5_qE;
elseif N==5
    %%%%%%%%%%EH at UAV%%%%%%%%%%%
    PU1 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    PU2 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    PU3 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    PU4 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
   % PU5 = tau*delta*PB*(N+1)*Omega_PB/(1-tau);
    %%%%%%%%%%%% OP phase 1 %%%%%%%%%%%%%%%%%%
    tu_1_P = exp(-Xi_O_P*d_GP^2*((PS+PP)*Omega_e+1)/(PP*Omega_GP));
    mau_1_P = Omega_SP*((Xi_O_P*d_GP^2*PS/(d_SP^2*PP*Omega_GP))+1/Omega_SP);
    O_1_P = 1-tu_1_P/mau_1_P;
    
    %%%%%%%%%%%% OP phase 2 %%%%%%%%%%%%%%%%%%
    tu_2_P = ((m/Omega_U1P)^m)*exp(-Xi_O_P*d_GP^2*((PU1+PP)*Omega_e+1)/(PP*Omega_GP))*factorial(m-1);
    mau_2_P = gamma(m)*((Xi_O_P*d_GP^2*PU1/(L_U1P*PP*Omega_GP))+m/Omega_U1P)^2;
    O_2_P = 1-tu_2_P/mau_2_P;
    
    %%%%%%%%%%%% OP phase n %%%%%%%%%%%%%%%%%%
    tu_3_P = exp(-Xi_O_P*d_GP^2*(PP*Omega_e+1)/(PP*Omega_GP));
    O_3_P = 1-tu_3_P;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P1_p = 0;
    P2_p = 0;
    P3_p = 0;
    P4_p = 0;
    P4_p_temp = 0;
   % P5_p_temp =0;
   % P5_p_tempE =0;
    for j = 0 : m-1
        A1 = m*Xi_Op*L_SU1*(PS*Omega_e + I_P +1)/(Omega_SU1*mu_p*PS);
        P1_p =P1_p + (A1^j)*exp(-A1)/factorial(j);
        
        A2 =  m*Xi_Op*L_U1U2*(PU1*Omega_e + I_P +1)/(Omega_U1U2*mu_p_1*PU1);
        P2_p =P2_p + (A2^j)*exp(-A2)/factorial(j);
        
        A3 =  m*Xi_Op*L_U2U3*(PU2*Omega_e + I_P +1)/(Omega_U2U3*mu_p_2*PU1);
        P3_p =P3_p + (A3^j)*exp(-A3)/factorial(j);
        
        A4 =  m*Xi_Op*L_U3U4*(PU3*Omega_e + I_P +1)/(Omega_U3U4*mu_p_3*PU1);
        P4_p =P4_p + (A4^j)*exp(-A4)/factorial(j);
        
        %%%%%%%%%%%%%%%
        Temp_3_p = 0;
        for k = 0:j
            A3_tu = nchoosek(j,k)*(((PP+PU4)*Omega_e + I_P +1)^(j-k))*((PP/L_GU5)^k)*factorial(k+m-1);
            A3_mau = gamma(m)*(((m*Xi_Op*L_U4U5*PP/(Omega_U4U5*mu_p_4*PU4*L_GU5))+m/Omega_GU5)^(k+m));
            Temp_3_p = Temp_3_p + A3_tu/A3_mau;
        end
        P4_p_temp = P4_p_temp + (Temp_3_p*(m*Xi_Op*L_U4U5/(Omega_U4U5*mu_p_4*PU4))^j)/factorial(j);
        %%%%%%%%%%%%
       % Temp_5_p = 0;
       % for k = 0:j
          %  A5_tu = nchoosek(j,k)*(((PP+PU5)*Omega_e + I_P +1)^(j-k))*((PP/d_GDp^2)^k)*factorial(k);
           % A5_mau = (m*Xi_Op*L_U5Dp*PP/(Omega_U5Dp**PU4*d_GDp^2)+1/Omega_GDp)^(k+1);
           % Temp_5_p = Temp_5_p+A5_tu/A5_mau;
       % end
        %P5_p_temp = P5_p_temp + (Temp_5_p*(m*Xi_Op*L_U5Dp/(Omega_U5Dp**PU5))^j)/factorial(j);
        %%%%%%%%%%%%EAV%%%%%%%%%%%%%
       % Temp_5_pE = 0;
       % for k = 0:j
         %   A5_tu = nchoosek(j,k)*(((PP+PU5)*Omega_e + I_P +1)^(j-k))*((PP/d_GE^2)^k)*factorial(k);
         %   A5_mau = (m*Xi_Op*L_U5E*PP/(Omega_U5E**PU4*d_GE^2)+1/Omega_GE)^(k+1);
         %   Temp_5_pE = Temp_5_pE+A5_tu/A5_mau;
       % end
        %P5_p_tempE = P5_p_tempE + (Temp_5_p*(m*Xi_E*L_U5E/(Omega_U5E**PU5))^j)/factorial(j);
        %%%%%%%%%%%%%%%%%%%%
    end
   % P5_p = P4_p_temp*((m/Omega_GU5)^m)*exp(-m*Xi_Op*L_U4U5*((PP+PU4)*Omega_e + I_P +1)/(Omega_U4U5*mu_p_4*PU4));
   % P6_p = P5_p_temp*exp(-(m*Xi_Op*L_U5Dp*((PP+PU5)*Omega_e + I_P +1))/(Omega_U5Dp**PU5))/Omega_GDp;
    
   % P6_pE = P5_p_tempE*exp(-(m*Xi_E*L_U5E*((PP+PU5)*Omega_e + I_P +1))/(Omega_U5E**PU5))/Omega_GE;
    
    O_Dp = 1-P1_p*P2_p*P3_p*P4_p*P5_p*P6_p;
    Th_Dp = (1-O_Dp)*R_Op;
    
    L_Dp = P1_p*P2_p*P3_p*P4_p*P5_p*P6_pE;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Throughput of the secondary network Dq
    P1_q = 0;
    P2_q = 0;
    P3_q = 0;
    P4_q = 0;
    P3_q_temp = 0;
    P4_q_temp =0;
    P4_q_tempE = 0;
    for j = 0 : m-1
        B1 = m*Xi_Oq*L_SU1*(PS*Omega_e+1)/(Omega_SU1*(mu_q-mu_p*Xi_Oq)*PS);
        P1_q =P1_q + (B1^j)*exp(-B1)/factorial(j);
        
        B2 =  m*Xi_Oq*L_U1U2*(PU1*Omega_e+1)/(Omega_U1U2*(mu_q_1-mu_p_1*Xi_Oq)*PU1);
        P2_q =P2_q + (B2^j)*exp(-B2)/factorial(j);
        
        B3 =  m*Xi_Oq*L_U2U3*(PU2*Omega_e+1)/(Omega_U2U3*(mu_q_2-mu_p_2*Xi_Oq)*PU2);
        P3_q =P3_q + (B3^j)*exp(-B3)/factorial(j);
        
        B4 =  m*Xi_Oq*L_U3U4*(PU3*Omega_e+1)/(Omega_U3U4*(mu_q_3-mu_p_3*Xi_Oq)*PU3);
        P4_q =P4_q + (B4^j)*exp(-B4)/factorial(j);
        
        %%%%%%%%%%%%%%%
        Temp_3_q = 0;
        for k = 0:j
            B3_tu = nchoosek(j,k)*(((PP+PU4)*Omega_e+1)^(j-k))*((PP/L_GU5)^k)*factorial(k+m-1);
            B3_mau = gamma(m)*(((m*Xi_Oq*L_U4U5*PP/(Omega_U4U5*(mu_q_4-mu_p_4*Xi_Oq)*PU4*L_GU5))+m/Omega_GU5)^(k+m));
            Temp_3_q = Temp_3_q + B3_tu/B3_mau;
        end
        P3_q_temp = P3_q_temp + (Temp_3_q*(m*Xi_Oq*L_U4U5/(Omega_U4U5*(mu_q_4-mu_p_4*Xi_Oq)*PU4))^j)/factorial(j);
        %%%%%%%%%%%%
        Temp_4_q = 0;
        for k = 0:j
            B4_tu = nchoosek(j,k)*(((PP+PU5)*Omega_e+1)^(j-k))*((PP/d_GDq^2)^k)*factorial(k);
            B4_mau = (m*Xi_Oq*L_U5Dq*PP/(Omega_U5Dq*(-*Xi_Oq)*PU5*d_GDq^2)+1/Omega_GDq)^(k+1);
            Temp_4_q = Temp_4_q+B4_tu/B4_mau;
        end
        P4_q_temp = P4_q_temp + (Temp_4_q*(m*Xi_Oq*L_U5Dq/(Omega_U5Dq*(-*Xi_Oq)*PU5))^j)/factorial(j);
        %%%%%%%%%%%%EAV%%%%%%%%%%%%%
        Temp_4_qE = 0;
        for k = 0:j
            B4_tu = nchoosek(j,k)*(((PP+PU5)*Omega_e+1)^(j-k))*((PP/d_GE^2)^k)*factorial(k);
            B4_mau = (m*Xi_E*L_U5E*PP/(Omega_U5E*(-*Xi_E)*PU5*d_GE^2)+1/Omega_GE)^(k+1);
            Temp_4_qE = Temp_4_qE+B4_tu/B4_mau;
        end
        P4_q_tempE = P4_q_tempE + (Temp_4_qE*(m*Xi_E*L_U5E/(Omega_U5E*(-*Xi_E)*PU5))^j)/factorial(j);
        %%%%%%%%%%%%%%%%%%%%%%
    end
    P5_q = P3_q_temp*exp(-(m*Xi_Oq*L_U4Dq*((PP+PU4)*Omega_e+1))/(Omega_U4Dq*(-*Xi_Oq)*PU4))/Omega_GDq;
    %P6_q = P4_q_temp*exp(-(m*Xi_Oq*L_U5Dq*((PP+PU5)*Omega_e+1))/(Omega_U5Dq*(-*Xi_Oq)*PU5))/Omega_GDq;
    
    P5_qE = P3_q_tempE*exp(-(m*Xi_E*L_U4E*((PP+PU4)*Omega_e+1))/(Omega_U34*(-*Xi_E)*PU4))/Omega_GE;
    
    O_Dq = 1-P1_q*P2_q*P3_q*P4_q*P5_q;
    Th_Dq = (1-O_Dq)*R_Oq;
    
    
    L_Dq = P1_q*P2_q*P3_q*P4_q*P5_qE;
    
% end
end

output = Th_Dp+Th_Dq;

function p_select = select(p,fitness)

for i = 1:length(fitness)
    h(i) = sum(fitness(1:i));
end

% y = y;
id_array = [];
for i = 1:length(fitness) %number of candidates
    maxx = h(length(h));
    r = rand(1)*maxx;
    comp = r > h;
    ids = find(comp==0);
    id = ids(1);
    id_array = [id_array id];
end
id_array;
p_select = p(:,id_array);


