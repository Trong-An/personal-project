function output = mophong(gamma_P,m,U1,h1,N,S,P,G,PP,Omega_e,R_O_P,W,nol)

d_GP = sqrt((G(1)-P(1))^2+(G(2)-P(2))^2);
d_SP = sqrt((S(1)-P(1))^2+(S(2)-P(2))^2);

d_U1P = sqrt((U1(1)-P(1))^2+(U1(2)-P(2)).^2+(h1).^2);
Xi_O_P = (2^((N+1)*R_O_P/W))-1;

% Best Relay Scheme
        
        gGP = max(gamrnd(1,1,1,1,nol),[],1);
        gSP = max(gamrnd(1,1,1,1,nol),[],1);
        gU1P = max(gamrnd(m,1/m,1,1,nol),[],1);
        
%        slot 1
        
        gamma_1_P = PP.*gGP./((d_GP^2).*((PP.*gSP./(d_SP^2))+(PP+PP)*Omega_e)+1);
        
%        slot 2
       
        gamma_2_P = PP.*gGP./((d_GP.^2).*((PP.*gU1P./(d_U1P.^2))+(PP+PP)*Omega_e)+1);
        
%        slot 3
        
        gamma_3_P = PP.*gGP./((d_GP.^2).*((PP.*gU1P./(d_U1P.^2))+(PP+PP)*Omega_e)+1);
        
  %      slot 4
        
        gamma_4_P = PP.*gGP./((d_GP^2).*(PP*Omega_e+1));
   %    Slot 5

        gamma_5_P = PP.*gGP./((d_GP^2).*(PP*Omega_e+1));
        
        gamma_P1 = min(gamma_1_P,gamma_2_P);
        gamma_P2 = min(gamma_3_P,gamma_4_P);
        gamma_P3 = min(gamma_P2,gamma_5_P);
        gamma_P = min(gamma_P1,gamma_P3);
   
countP = 0;
        
            %  for j = 1:length(PPdB)
            for  i = 1:nol
                if gamma_P(i) <= Xi_O_P 
                    countP = countP +1;
                    
                   output = countP/nol;
                   
                    %countP
                   % O_P
                end
            end
            %  end
           % countP
           % O_P
            % output = O_P;
