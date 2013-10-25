clc;
SIZE=1001;
maxTime = 5001;

%Constants
c=299792458;
% Medium Specifications
mu=[1.2566e-006*ones(1,SIZE)];   %permeability of free sapce
% epsilon=[8.8542e-012*ones(1,SIZE)]; %permittivity of free space 
epsilon=[8.8542e-012*ones(1,SIZE-180) 1.7708e-011*ones(1,120) 8.8542e-012*ones(1,160)]; %introducing medium of 2*e with a width of 20

% Courant Number (Accuracy) Sc
% ideal Condition --> Sc= c*delt/delx = 1
% f=3Ghz, lambda=c/f=0.1m, for 4 wavelengths, dx=0.4/(maxtime=1000)
f=3e9;
lambda=c/f;
delx=(4*lambda)/SIZE;
% so dt=dx/c=1.333e-12
delt=delx/c;
Sc=c*delt/delx;
epsilonr=1;
mur=1;

% Temp Variable
ez=zeros(1,SIZE);
hy=zeros(1,SIZE-1);
mm=0;
temp=0;
ez1q=0;
ez2q=0;
ezmq=0;
ezm1q=0;

    for qTime = 1:(maxTime-1)
%        Update Magnetic field
        for  mm = 1:(SIZE-1)
            hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) * (delt/(delx*mu(mm)));
        end
%         Update Electrical filed
        for mm = 2:(SIZE-1)
            ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * (delt/(delx*epsilon(mm))) ;
        end
%         Source node (hard coded)
        ez(2) = ez(2)+ (sin(2*pi*(qTime)*f*delt)*Sc);
%         Absorbing Boundary Conditions
        ez(1)=ez2q+(ez(2)-ez1q)*(((Sc/(mur*(epsilonr))^0.5)-1)/((Sc/(mur*(epsilonr))^0.5)+1));
        ez(SIZE)=ezm1q+(ez(SIZE-1)-ezmq)*(((Sc/(mur*(epsilonr))^0.5)-1)/((Sc/(mur*(epsilonr))^0.5)+1));
%         Saving q-1 (pervious step time values) for boundary Conditions
        ez2q=ez(2);
        ez1q=ez(1);
        ezmq=ez(SIZE);
        ezm1q=ez(SIZE-1);
%         Plotting
        figure(1);
        grid on; 
        subplot(2,1,1);
        plot(1:SIZE,ez);
        title('Electirc Component');
        xlim([0 SIZE]);
        ylim([-1.2 1.2]);
        line([SIZE-180 SIZE-180],[-1.2 1.2],'Color','Red') % Medium slab line
        line([SIZE-160 SIZE-160],[-1.2 1.2],'Color','Red') % Medium slab line
        subplot(2,1,2);
        plot(1:SIZE-1,hy);
        title('Magnetic Component');
        xlim([0 SIZE]);
        ylim([-0.005 0.005]);
        line([SIZE-180 SIZE-180],[-0.005 0.005],'Color','Red') % Medium slab line
        line([SIZE-160 SIZE-160],[-0.005 0.005],'Color','Red') % Medium slab line
%         pause(0.02);
    end
