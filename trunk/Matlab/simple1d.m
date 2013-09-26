SIZE=201;
ez=zeros(1,SIZE);
hy=zeros(1,SIZE);
imp0=377; %squareroot(u0/e0)
maxTime = 1001;
mm=0;
temp=0;
    for qTime = 1:(maxTime-1)
%        Update Magnetic field
        for  mm = 1:(SIZE - 2)
            hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) / imp0;
        end
%         Update Electrical filed
        for mm = 2:(SIZE-1)
            ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * imp0;
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
        ylim([-1.2 1.2]);
        subplot(2,1,2);
        plot(1:SIZE,hy);
        title('Magnetic Component');
        Xlim([0 SIZE]);
        ylim([-0.003 0.003]);
        pause(0.02);
    end
