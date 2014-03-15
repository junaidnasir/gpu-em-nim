SIZE=201;
ez=zeros(1,SIZE);
hy=zeros(1,SIZE);
mu=[1.2566e-006*ones(1,SIZE)];   %permeability of free sapce
% epsilon=[8.8542e-012*ones(1,SIZE)]; %permittivity of free space 
epsilon=[8.8542e-012*ones(1,SIZE-80) 1.7708e-011*ones(1,20) 8.8542e-012*ones(1,60)]; %introducing medium of 2*e with a width of 20
% ideal Condition --> Sc= c*delt/delx = 1
delt=1;
delx=299792458;
maxTime = 1001;
mm=0;
temp=0;
    for qTime = 1:(maxTime-1)
%        Update Magnetic field
        for  mm = 1:(SIZE - 2)
            hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) * (delt/(delx*mu(mm)));
        end
%         Update Electrical filed
        for mm = 2:(SIZE-1)
            ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * (delt/(delx*epsilon(mm))) ;
        end
%         Source node (hard coded)
        ez(1) = exp(-(qTime - 30) * (qTime - 30) / 100.);
%         Plotting
        figure(1);
        grid on; 
        subplot(2,1,1);
        plot(1:SIZE,ez);
        title('Electirc Component');
        Xlim([0 SIZE]);
        ylim([-1.5 1.5]);
        line([SIZE-80 SIZE-80],[-1.5 1.5],'Color','Red')
        line([SIZE-60 SIZE-60],[-1.5 1.5],'Color','Red')
        subplot(2,1,2);
        plot(1:SIZE,hy);
        title('Magnetic Component');
        Xlim([0 SIZE]);
        ylim([-0.005 0.005]);
        line([SIZE-80 SIZE-80],[-0.005 0.005],'Color','Red')
        line([SIZE-60 SIZE-60],[-0.005 0.005],'Color','Red')
        pause(0.02);
    end
