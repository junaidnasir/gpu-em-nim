SIZE=201;
ez=zeros(1,SIZE);
hy=zeros(1,SIZE);
imp0=377;
maxTime = 1001;
mm=0;
temp=0;
    for qTime = 1:(maxTime-1)

        for  mm = 1:(SIZE - 2)
            hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) / imp0;
        end
        for mm = 2:(SIZE-1)
            ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * imp0;
        end
        ez(1) = exp(-(qTime - 30) * (qTime - 30) / 100.);
        figure(1);
        grid on; 
        ylim([-1.2 1.2])
        plot(1:201,ez);
        pause(0.02);
    end


