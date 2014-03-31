tic
clc;
SIZE=5120;
maxTime = 1024;
SourceSelect=1; % 0=Sinosoidal, 1=Gauassian
if (SourceSelect==0)
    maxTime = 5001;
end
%Constants
c=3e8;

% Courant Number (Accuracy) Sc
% ideal Condition --> Sc= c*delt/delx = 1
% f=3Ghz, lambda=c/f=0.1m, for 4 wavelengths, dx=0.4/(maxtime=1000)
PulseWidth=800;
f=3e9;
w=2*pi*f;    % omega
k0=w/c     ; % free space wave number constant
lambda=c/f;
delx=(4*lambda)/SIZE;
% so dt=dx/c=1.333e-12
delt=delx/c
Sc=c*delt/delx;
epsilonr=1;
mur=1;


% Incident and Refelected Waves Variables
Eincident=0;
Etransmitted=0;
Etemp=zeros(1,maxTime);


% refractive index variables
Z1 = 750;
z1 = Z1*delx;
Z2 = 760;
z2 = Z2*delx;
Exz1 = zeros(maxTime,1); % record Electric field at 750
Exz2 = zeros(maxTime,1); % record electric field at 760
fspan = 100; % Points to plot in frequency domain

for medium= 1:2
    % Temp Variable
    ez=zeros(1,SIZE);
    hy=zeros(1,SIZE-1);
    mm=0;
    ez1q=0;
    ez2q=0;
    ezmq=0;
    ezm1q=0;
    % Medium Specifications
    mu=1.2566e-006*ones(1,SIZE);   %permeability of free sapce
    if medium==1
        epsilon=8.8542e-012*ones(1,SIZE); % free space
    else
        epsilon=[8.8542e-012*ones(1,SIZE-500) 1.7708e-011*ones(1,500)]; % half medium
    end
    for qTime = 1:(maxTime-1)
%        Update Magnetic field
    	for  mm = 1:(SIZE-1)
            hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) * (delt/(delx*mu(mm)));
        end
%         Update Electrical filed
        for mm = 2:(SIZE-1)
            ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * (delt/(delx*epsilon(mm))) ;
        end
        Etemp(qTime)= ez(SIZE-498); %just after boundary of medium
        if SourceSelect==0
%         Source node (hard coded)
		    ez(2) = ez(2)+ (sin(2*pi*(qTime)*f*delt)*Sc);
		else
		    ez(2) = ez(2)+exp(-(qTime - 30) * (qTime - 30) / (PulseWidth./4));
		end
%         Absorbing Boundary Conditions
        ez(1)=ez2q+(ez(2)-ez1q)*(((Sc/(mur*(epsilonr))^0.5)-1)/((Sc/(mur*(epsilonr))^0.5)+1));
        ez(SIZE)=ezm1q+(ez(SIZE-1)-ezmq)*(((Sc/(mur*(epsilonr))^0.5)-1)/((Sc/(mur*(epsilonr))^0.5)+1));
%         Saving q-1 (pervious step time values) for boundary Conditions
		ez2q=ez(2);
		ez1q=ez(1);
		ezmq=ez(SIZE);
		ezm1q=ez(SIZE-1);
%         Plotting
%         figure(1);
%         subplot(2,1,1);
%         plot(1:SIZE,ez);
%         title('Electirc Component');
%         xlim([0 SIZE]);
%         ylim([-1.2 1.2]);
%         if medium==2
%             line([SIZE-500 SIZE-500],[-1.2 1.2],'Color','Red') % Medium slab line
%         end
%         subplot(2,1,2);
%         plot(1:SIZE-1,hy);
%         title('Magnetic Component');
%         xlim([0 SIZE]);
%         ylim([-0.005 0.005]);
%         if medium==2
%         line([SIZE-500 SIZE-500],[-0.005 0.005],'Color','Red') % Medium slab line
%         end
          Exz1(qTime)=ez(Z1);
          Exz2(qTime)=ez(Z2);
    end
    if medium==1
        Eincident=Etemp;
    else
        Etransmitted=Etemp;
    end
end
toc
% % Frequency Domain Analysis
% Fs=1/delt;   %Sampling Frequency
% T=1/Fs;      %Sample Time
% L=maxTime;      %Length of Signal
% time=(0:L-1)*T;    %Time Vector
% figure(2);
% subplot(2,1,1);
% plot(Fs*time,Eincident)
% title('Incident Wave');
% xlabel('time (picoseconds)')
% subplot(2,1,2);
% plot(Fs*time,Etransmitted);
% title('Transmitted Wave');
% xlabel('time (picoseconds)')
% % Fourrier Domain
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% FEincident = fft(Eincident,NFFT)/L;
% FEtransmitted = fft(Etransmitted,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);          %frequency scaling
% % Plot single-sided amplitude spectrum.
% figure(3);
% subplot(2,1,1);
% plot(f,2*abs(FEincident(1:NFFT/2+1))) 
% xlim([0 1e11]);
% % ylim([0 0.5]);
% title('Single-Sided Amplitude Spectrum of Incident Wave')
% xlabel('Frequency (Hz)')
% ylabel('|Eincident(f)|')
% subplot(2,1,2);
% plot(f,2*abs(FEtransmitted(1:NFFT/2+1)))
% xlim([0 1e11]);
% % ylim([0 0.5]);
% title('Single-Sided Amplitude Spectrum of Transmitted Wave')
% xlabel('Frequency (Hz)')
% ylabel('|Etransmitted(f)|')
% 
% cTransmitted=FEtransmitted/FEincident
% cReflected=1-cTransmitted
% eta1=sqrt(1/1);
% eta2=sqrt(2/1);
% Gamma=(eta2-eta1)/(eta2+eta1)
% 
% 
% % eq 33
% EXZ1 = fft(Exz1,NFFT)/L;
% EXZ2 = fft(Exz2,NFFT)/L;
% 
% nFDTD = (1/(1i*k0*(z1-z2))).*log(EXZ2(1:NFFT/2+1)./EXZ1(1:NFFT/2+1));
% figure(4);
% plot(f(1:fspan), real(nFDTD(1:fspan)));
% title('Refractive index re(n)');
% xlabel('Frequency (Hz)');
% ylabel('re(n)');
% line([3e9 3e9],[-15 1.415],'Color','Red')
% line([0e9 3e9],[1.415 1.415],'Color','Red')
% 
% ReferectiveIndex=(1/(k0*(760-750)*i))*log(FEtransmitted(760)/FEtransmitted(750))