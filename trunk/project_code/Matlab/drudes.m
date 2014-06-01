%% Variables and stuff
clear;
clc;
a=13;
SIZE=1048;   
maxTime =5*SIZE;

SourceSelect=0; % 0=Sinosoidal, 1=Gauassian
% if (SourceSelect==0)
%     maxTime = 5001;
% end
%Constants
c=3e8;

% Courant Number (Accuracy) Sc
% ideal Condition --> Sc= c*delt/delx = 1
% f=3Ghz, lambda=c/f=0.1m, for 4 wavelengths, dx=0.4/(maxtime=1000)
PulseWidth=800;
f=3e9;      % 3GHz
w=2*pi*f;    % omega
k0=w/c     ; % free space wave number constant
lambda=c/f;
delx=(4*lambda)/SIZE;
% so dt=dx/c=1.333e-12
delt=delx/c;
Sc=c*delt/delx;
epsilonr=-1;
mur=-1;
mu_nought=1.2566e-006;
epsilon_nought=8.8542e-012;
% Drudes model variables
mu_inf=1;
    epsilon_inf=1;
    omega_pe2= 4*pi*pi*f*f*(epsilon_inf-epsilonr);%Plasma frequency electric squared
    omega_pm2= 4*pi*pi*f*f*(mu_inf-mur);%Plasma frequency magnetic squared
    ro_m=0;
    ro_e=0;
% Incident and Refelected Waves Variables
Eincident=0;
Etransmitted=0;
Etemp=zeros(1,maxTime);


% refractive index variables
Z1 = (SIZE/2)+100;
z1 = Z1*delx;
Z2 = (SIZE/2)+110;
z2 = Z2*delx;
Exz1 = zeros(maxTime,1); % record Electric field at 750
Exz2 = zeros(maxTime,1); % record electric field at 760
fspan = 100; % Points to plot in frequency domain
%% Main Loop
for medium= 1:2
    % Medium Specifications
    if medium==1
        mu=1.2566e-006*ones(1,SIZE);   %permeability of free sapce
        epsilon=8.8542e-012*ones(1,SIZE); % free space
        epsilonr=ones(1,SIZE);
        mur=ones(1,SIZE);
    else
        epsilon=[8.8542e-012*ones(1,SIZE-(SIZE/2)) -8.8542e-012*ones(1,(SIZE/4)) 8.8542e-012*ones(1,(SIZE/4))]; % half medium
        mu=[1.2566e-006*ones(1,SIZE-(SIZE/2)) -1.2566e-006*ones(1,(SIZE/4)) 1.2566e-006*ones(1,(SIZE/4))];
        epsilonr=[ones(1,SIZE-(SIZE/2)) -1*ones(1,SIZE/4) ones(1,SIZE/4)];
        mur=[ones(1,SIZE-(SIZE/2)) -1*ones(1,SIZE/4) ones(1,SIZE/4)];
    end
    for  mm = 1:(SIZE)
        omega_pe2= 4*pi*pi*f*f*(epsilon_inf-epsilonr(mm));%Plasma frequency electric squared
        omega_pm2= 4*pi*pi*f*f*(mu_inf-mur(mm));%Plasma frequency magnetic squared
        m_divide=((4*mu_nought*mu_inf)+(mu_nought*omega_pm2*(delt*delt))+(mu_nought*mu_inf*ro_m*(2*delt)));
        am(mm)= 4/m_divide;
        bm(mm)= (ro_m*2*delt)/m_divide;
        cm(mm)= (4*mu_nought*mu_inf)/m_divide;
        dm(mm)= (-mu_nought*omega_pm2*delt*delt)/m_divide;
        em(mm)= (mu_nought*mu_inf*ro_m*2*delt)/m_divide;
        e_divide=((4*epsilon_nought*epsilon_inf)+(epsilon_nought*omega_pe2*(delt*delt))+(epsilon_nought*epsilon_inf*ro_e*(2*delt)));
        ae(mm)= 4/e_divide;
        be(mm)= (ro_e*2*delt)/e_divide;
        ce(mm)= (4*epsilon_nought*epsilon_inf)/e_divide;
        de(mm)= (-epsilon_nought*omega_pe2*delt*delt)/e_divide;
        ee(mm)= (epsilon_nought*epsilon_inf*ro_e*2*delt)/e_divide;
    end
    % Temp Variable
    ez=zeros(1,SIZE);
    ezn_0=zeros(1,SIZE-1);
    ezn_1=zeros(1,SIZE-1);
    hy=zeros(1,SIZE-1);
    hyn_0=zeros(1,SIZE-1);
    hyn_1=zeros(1,SIZE-1);
    by=zeros(1,SIZE);
    byn_0=zeros(1,SIZE);
    byn_1=zeros(1,SIZE);
    dz=zeros(1,SIZE);
    dzn_0=zeros(1,SIZE);
    dzn_1=zeros(1,SIZE);
    mm=0;
    ez1q=0;
    ez2q=0;
    ezmq=0;
    ezm1q=0;
    
    for qTime = 1:(maxTime-1)
        % Drudes model Calculations
        
        
        %storing pervious time steps values
        hyn_1=hyn_0;
        hyn_0=hy;
        byn_1=byn_0;
        byn_0=by;
        ezn_1=ezn_0;
        ezn_0=ez;
        dzn_1=dzn_0;
        dzn_0=dz;
%        Update By
        for  mm = 1:(SIZE-1)
            by(mm) = by(mm) + ((ez(mm) - ez(mm+1)) * (delt/(delx)));
        end
%        Update Magnetic field
        for  mm = 1:(SIZE-1)  %changed it from 1 to 2 because of by(0) index problem at 1
            hy(mm) = (am(mm)*(by(mm)-2*byn_0(mm)+byn_1(mm)))+(bm(mm)*(by(mm)-byn_1(mm)))+(cm(mm)*((2*hyn_0(mm))-(hyn_1(mm))))+(dm(mm)*((2*hyn_0(mm))+(hyn_1(mm))))+(em(mm)*(hyn_1(mm)));
        end
%         Update Dz
        for mm = 2:(SIZE-1)
            dz(mm) = dz(mm) + ((hy(mm-1) - hy(mm)) * (delt/delx));
        end
%         Update Electrical filed
        for mm = 2:(SIZE-1)
            ez(mm) = (ae(mm)*(dz(mm)-(2*dzn_0(mm))+dzn_1(mm)))+(be(mm)*(dz(mm)-dzn_1(mm)))+(ce(mm)*((2*ezn_0(mm))-ezn_1(mm)))+(de(mm)*((2*ezn_0(mm))+ezn_1(mm)))+(ee(mm)*ezn_1(mm));
        end
        Etemp(qTime)= ez(SIZE-(SIZE/2)+2); %just after boundary of medium
        if SourceSelect==0
%         Source node (hard coded)
		    ez(2) = ez(2)+ (sin(2*pi*(qTime)*f*delt)*Sc);
         else
		    ez(2) = ez(2)+exp(-(qTime - 30) * (qTime - 30) / (PulseWidth./4));
        end
%         Absorbing Boundary Conditions
        ez(1)=ez2q+(ez(2)-ez1q)*(((Sc/(mur(1)*(epsilonr(1)))^0.5)-1)/((Sc/(mur(1)*(epsilonr(1)))^0.5)+1));
        ez(SIZE)=ezm1q+(ez(SIZE-1)-ezmq)*(((Sc/(mur(1)*(epsilonr(1)))^0.5)-1)/((Sc/(mur(1)*(epsilonr(1)))^0.5)+1));
        dz(SIZE)=epsilonr(1)*epsilon_nought; %epsilon nough
%         Saving q-1 (pervious step time values) for boundary Conditions
		ez2q=ez(2);
		ez1q=ez(1);
		ezmq=ez(SIZE);
		ezm1q=ez(SIZE-1);
%         Plotting
       if mod(qTime,10)==0
       if medium==2
        figure(1);
%         subplot(2,1,1);
        plot(1:SIZE,ez);
        title('Electirc Component');
        xlim([0 SIZE]);
%         ylim([-1.2 1.2]);
%         if medium==2
            line([SIZE-(SIZE/2) SIZE-(SIZE/2)],[-50 50],'Color','Red') % Medium slab line
            line([SIZE-(SIZE/4) SIZE-(SIZE/4)],[-50 50],'Color','Red') % Medium slab line
%         end
       end
       end
          Exz1(qTime)=ez(Z1);
          Exz2(qTime)=ez(Z2);
          
       
    end
    if medium==1
        Eincident=Etemp;
    else
        Etransmitted=Etemp;
    end
end
%% Post processing
% Frequency Domain Analysis
Fs=1/delt;   %Sampling Frequency
T=1/Fs;      %Sample Time
L=maxTime;      %Length of Signal
time=(0:L-1)*T;    %Time Vector
figure(2);
subplot(2,1,1);
plot(Fs*time,Eincident)
title('Incident Wave');
xlabel('time (picoseconds)')
subplot(2,1,2);
plot(Fs*time,Etransmitted);
title('Transmitted Wave');
xlabel('time (picoseconds)')
% Fourrier Domain
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
FEincident = fft(Eincident,NFFT)/L;
FEtransmitted = fft(Etransmitted,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);          %frequency scaling
% Plot single-sided amplitude spectrum.
figure(3);
subplot(2,1,1);
plot(f,2*abs(FEincident(1:NFFT/2+1))) 
xlim([0 1e11]);
% ylim([0 0.5]);
title('Single-Sided Amplitude Spectrum of Incident Wave')
xlabel('Frequency (Hz)')
ylabel('|Eincident(f)|')
subplot(2,1,2);
plot(f,2*abs(FEtransmitted(1:NFFT/2+1)))
xlim([0 1e11]);
% ylim([0 0.5]);
title('Single-Sided Amplitude Spectrum of Transmitted Wave')
xlabel('Frequency (Hz)')
ylabel('|Etransmitted(f)|')

cTransmitted=FEtransmitted/FEincident
cReflected=1-cTransmitted
eta1=sqrt(1/1);
eta2=sqrt(2/1);
Gamma=(eta2-eta1)/(eta2+eta1)


% eq 33
EXZ1 = fft(Exz1,NFFT)/L;
EXZ2 = fft(Exz2,NFFT)/L;

nFDTD = (1/(1i*k0*(z1-z2))).*log(EXZ2(1:NFFT/2+1)./EXZ1(1:NFFT/2+1));
figure(4);
plot(f(1:fspan), real(nFDTD(1:fspan)));
title('Refractive index re(n)');
xlabel('Frequency (Hz)');
ylabel('re(n)');
line([3e9 3e9],[-15 -1],'Color','Red')
line([0e9 3e9],[-1 -1],'Color','Red')

ReferectiveIndex=(1/(k0*(Z2-Z1)*i))*log(FEtransmitted(Z2)/FEtransmitted(Z1))