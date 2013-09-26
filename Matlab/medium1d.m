SIZE=201;
ez=zeros(1,SIZE);
hy=zeros(1,SIZE);
% introducing a medium having e=2*e0 (slab of length 20) -->
% imp=squareroot(u0/2*e0)=265.77
imp=[377*ones(1,SIZE-80) 266*ones(1,20) 377*ones(1,60)];
% imp0=377; %squareroot(u0/e0) = 377
maxTime = 1001;
mm=0;
temp=0;
    for qTime = 1:(maxTime-1)
%        Update Magnetic field
        for  mm = 1:(SIZE - 2)
            hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) / imp(1); %stable
%             with some refelction of waves
%             hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) / imp(mm); %unstable
        end
%         Update Electrical filed
        for mm = 2:(SIZE-1)
            ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * imp(mm);
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
%         ylim([-1.2 1.2]);
        z=ylim;
        line([SIZE-80 SIZE-80],[z(1) z(2)],'Color','Red')
        line([SIZE-60 SIZE-60],[z(1) z(2)],'Color','Red')
        subplot(2,1,2);
        plot(1:SIZE,hy);
        title('Magnetic Component');
        Xlim([0 SIZE]);
%         ylim([-0.003 0.003]);
        z=ylim;
        line([SIZE-80 SIZE-80],[z(1) z(2)],'Color','Red')
        line([SIZE-60 SIZE-60],[z(1) z(2)],'Color','Red')
        pause(0.02);
    end
