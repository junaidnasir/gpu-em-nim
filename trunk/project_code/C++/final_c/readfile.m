SIZE=1001;
medium=1;
fspan=100;
filepath=fullfile(pwd, 'results');

fidp = fopen(strcat(filepath,'\maxTime','.jd'),'r','l');
if fidp==-1
    display('Error');
    return
end
maxTime=fread(fidp,1,'int');
fclose(fidp);

fidp = fopen(strcat(filepath,'\data','.jd'),'r','l');
if fidp==-1
    display('Error');
    return
end
data=fread(fidp,4,'double');
fclose(fidp);
delt=data(1);
k0=data(2);
z1=data(3);
z2=data(4);

fidp = fopen(strcat(filepath,'\Eincident','.jd'),'r','l');
if fidp==-1
    display('Error');
    return
end
Eincident=fread(fidp,maxTime,'double');
fclose(fidp);

fidp = fopen(strcat(filepath,'\Etransmitted','.jd'),'r','l');
if fidp==-1
    display('Error');
    return
end
Etransmitted=fread(fidp,maxTime,'double');
fclose(fidp);

% Frequency Domain Analysis
Fs=1/delt;   %Sampling Frequency
T=1/Fs;      %Sample Time
L=maxTime;      %Length of Signal
time=(0:L-1)*T;    %Time Vector
figure(1);
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
figure(2);
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
 
cTransmitted=FEtransmitted/FEincident;
cReflected=1-cTransmitted;
eta1=sqrt(1/1);
eta2=sqrt(2/1);
Gamma=(eta2-eta1)/(eta2+eta1) 
 
% eq 33
fidp = fopen(strcat(filepath,'\Exz1','.jd'),'r','l');
if fidp==-1
    display('Error');
    return
end
Exz1=fread(fidp,maxTime,'double');
fclose(fidp);

fidp = fopen(strcat(filepath,'\Exz2','.jd'),'r','l');
if fidp==-1
    display('Error');
    return
end
Exz2=fread(fidp,maxTime,'double');
fclose(fidp);

EXZ1 = fft(Exz1,NFFT)/L;
EXZ2 = fft(Exz2,NFFT)/L;

nFDTD = (1/(1i*k0*(z1-z2))).*log(EXZ2(1:NFFT/2+1)./EXZ1(1:NFFT/2+1));
figure(3);
plot(f(1:fspan), real(nFDTD(1:fspan)));
title('Refractive index re(n)');
xlabel('Frequency (Hz)');
ylabel('re(n)');
line([3e9 3e9],[-15 1.415],'Color','Red')
line([0e9 3e9],[1.415 1.415],'Color','Red')

for medium= 1:2
    j=0;
    for qTime = 1:(maxTime-1)
        fidp = fopen(strcat(filepath,'\Efield',int2str(medium),'_',int2str(j),'.jd'),'r','l');
        if fidp==-1
            display('Error');
            return
        end
        ez=fread(fidp,maxTime,'double');
        fclose(fidp);
        j=j+1;
% Plotting
        figure(4);
        subplot(2,1,1);
        plot(1:SIZE,ez);
        title('Electirc Component');
        xlim([0 SIZE]);
        ylim([-1.2 1.2]);
        if medium==2
            line([SIZE-500 SIZE-500],[-1.2 1.2],'Color','Red') % Medium slab line
        end
        subplot(2,1,2);
%         plot(1:SIZE-1,hy);
        title('Magnetic Component');
        xlim([0 SIZE]);
        ylim([-0.005 0.005]);
        if medium==2
        line([SIZE-500 SIZE-500],[-0.005 0.005],'Color','Red') % Medium slab line
        end
    end
end