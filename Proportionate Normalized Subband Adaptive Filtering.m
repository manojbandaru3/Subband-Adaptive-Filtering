function [err] = PNSAF(N,Nw,dn,mu)

Nr = 100;
M  = 4;
u  = mu;
gamma = dn;
hk = FIR(M);
B  = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
A  = 1;
L  = M;
rho= 0.05;

for l=1:Nr
    disp(['Run: ' num2str(l)])
    nc=normrnd(0,0.01,1,N*M);       % noise at system output 
    xin = color(N*M,1);             % input signal
    d=filter(B,A,xin(1:end));       % desired signal
    d=d+nc;                         % unknown system output
    
    
     % Analysis of desired and input signals
    for k=1:M
        x_temp=filter(hk(k,:),1,d);
        dsb(k,:)=x_temp(find(mod((1:length(x_temp))-1,L)==0)); % desired signal split in subbands
        x_temp=filter(hk(k,:),1,xin);
        xsb(k,:)=x_temp(find(mod((1:length(x_temp))-1,L)==0)); % input signal split in subbands
    end;
    
    %initializing temporary values
    w_ol =zeros(1,Nw+1);      % initial coefficient vector for subband m
    for m=1:M
        x_ol(m,:)=[zeros(1,Nw+1)];    % initial input vector for subband m
    end
    
     %Adaption Process
    for k=1:N
      if mod(k,200)==0 disp(['Iteration: ' num2str(k)]), end;
      temp_error = zeros(1,Nw+1);
      for m=1:M 
          
          x_ol(m,:)=[xsb(m,k) x_ol(m,1:Nw)]; % new input vector for subband m
          x_temp=shiftdim(x_ol(m,:),1);
          w_temp=shiftdim(w_ol,1);
          e_ol(m,l,k)=dsb(m,k)-w_temp'*x_temp;   % error at subband m
          t = [];
          for j = 1:(Nw+1)
            b = [];
            b = [0.01,abs(w_ol)];
            tem = max(rho*max(b),(abs(w_ol(j))));
            t = [t,tem];
          end
          g = zeros(1,(Nw+1));
          for f = 1:(Nw+1)
              g(f) = ((Nw+1)*t(f)/sum(t));
          end
          G = diag(g);
          denom=((gamma/(Nw+1))+x_ol(m,:)*G*x_temp);
          temp_error = temp_error+( (e_ol(m,l,k)*shiftdim(G*x_temp,1))/denom);   
      end
      w_ol = w_ol + u*(temp_error); 
      
    end

end
for m=1:M
    mse_ol(m,:)=mean(e_ol(m,:,:).^2,2);    % MSE at subband m
end;
MSE=mean(mse_ol,1);                 % overall MSE
MSEsub = mse_ol;

err = 10*log10(MSE);
end