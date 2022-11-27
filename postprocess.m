clc; close all;

fprintf('Shape of Fin:\t\t\t\t %s\n', shape)
fprintf('Lamda Distribution:\t\t\t %s\n', heat_conduc)
fprintf('Simulation Type:\t\t\t %s\n', simulationType)

switch simulationType
    
    case 'steady'

        %% 2D Plot
        figure(1)
        pcolor(X,Y,T)
        hold on
        pcolor(X,-Y,T)
        if strcmp(shape, 'linear')
            axis([0 l -h1/2 h1/2 -1 1]) %limits of x & y axes
        elseif strcmp(shape, 'quadratic')
            axis([0 l -h1/2*1.15 h1/2*1.15 -1 1]) %limits of x & y axes
        else
            axis([0 l -h1/2*1.45 h1/2*1.45 -1 1]) %limits of x & y axes
        end
        set(gcf, 'Position',[70,150,550,550]);

        %% Contour Plot
        figure(2)
        contour(X,Y,T)
        hold on
        contour(X,-Y,T)
        set(gcf, 'Position',[635,150,550,550]);

        %% 3D Plot
        figure(3)
        surf(X,Y,T); hold on; surf(X,-Y,T)
        set(gcf, 'Position',[1200,150,550,550]);
        view(45,24);            %view angle


    case 'unsteady'

        if strcmp(TimeIntegrType, 'Theta') && theta == 0.5
            fprintf('Time Discretization Scheme:\t Crank-Nicolson\n')
        elseif strcmp(TimeIntegrType, 'Theta') && theta == 0
            fprintf('Time Discretization Scheme:\t Explicit\n')
        elseif strcmp(TimeIntegrType, 'Theta') && theta == 1
            fprintf('Time Discretization Scheme:\t Implicit\n')
        else
            fprintf('Time Discretization Scheme:\t %s\n', TimeIntegrType)
        end

        i=1; % plotting step
        if strcmp(TimeIntegrType, 'Explicit') || strcmp(TimeIntegrType, 'RungeKutta4')
            i=10;  % raise the step because those 2 methods need too small Δt
        end

        %% 3D Plot
        filename = 'Transient_3D.gif';
        figure(1)
        for t=1:i:size(T,3)
            surf(X,Y,T(:,:,t))
            hold on
            surf(X,-Y,T(:,:,t))
            hold off
            set(gcf, 'Position',[1200,150,550,550]);
            view(45,24);            %view angle

            % Record into a GIF
            drawnow
            frame= getframe(1);
            im= frame2im(frame);
            [imind,cm] = rgb2ind(im,64);
            if t == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
            end
        end

        %% 2D Plot
        filename = 'Transient_2D.gif';
        figure(2)
        for t=1:i:size(T,3)
            pcolor(X,Y,T(:,:,t))
            hold on
            pcolor(X,-Y,T(:,:,t))
            hold off
            if strcmp(shape, 'linear')
                axis([0 l -h1/2 h1/2 -1 1]) %limits of x & y axes
            elseif strcmp(shape, 'quadratic')
                axis([0 l -h1/2*1.15 h1/2*1.15 -1 1]) %limits of x & y axes
            else
                axis([0 l -h1/2*1.45 h1/2*1.45 -1 1]) %limits of x & y axes
            end
            set(gcf, 'Position',[70,150,550,550]);

            %Record into a GIF
            drawnow
            frame= getframe(2);
            im= frame2im(frame);
            [imind,cm] = rgb2ind(im,64);
            if t == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
            end
        end

        %% Contour plot
        filename = 'Transient_Contour.gif';
        figure(3)
        for t=1:i:size(T,3)
            contour(X,Y,T(:,:,t))
            hold on
            contour(X,-Y,T(:,:,t))
            hold off
            set(gcf, 'Position',[635,150,550,550]);

            %Record into a GIF
            drawnow
            frame= getframe(3);
            im= frame2im(frame);
            [imind,cm] = rgb2ind(im,64);
            if t == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
            end
        end

end