clc; close all;

fprintf('Shape of Fin:\t\t\t\t %s\n', shape)
fprintf('Lamda Distribution:\t\t\t %s\n', heat_conduc)
fprintf('Simulation Type:\t\t\t %s\n', simulationType)

switch simulationType
    
    case 'steady'

        %% 2D Plot
        figure(1)
        s = pcolor(X,Y,T);        hold on;
        s.FaceColor = 'interp';
        s = pcolor(X,-Y,T);       hold off;
        s.FaceColor = 'interp';
        if strcmp(shape, 'linear')
            axis([0 l -h1/2 h1/2 -1 1]) %limits of x & y axes
        elseif strcmp(shape, 'quadratic')
            axis([0 l -h1/2*1.15 h1/2*1.15 -1 1]) %limits of x & y axes
        else
            axis([0 l -h1/2*1.45 h1/2*1.45 -1 1]) %limits of x & y axes
        end
        xlabel('x'); ylabel('y');
        hcb=colorbar;
        set(get(hcb,'Title'),'String','T');     % insert label above colorbar
        set(gcf, 'Position',[10,150,620,550]);
        saveas(gcf,'Steady_2D.jpg')

        %% Contour Plot
        figure(2)
        contour(X,Y,T);        hold on;
        contour(X,-Y,T);       hold on;
        plot(X(1,:),Y(1,:),'k');    hold on;
        plot(X(1,:),-Y(1,:),'k');   hold off;
        xlabel('x'); ylabel('y');
        hcb=colorbar;
        set(get(hcb,'Title'),'String','T');
        set(gcf, 'Position',[640,150,620,550]);
        saveas(gcf,'Steady_Contour.jpg')

        %% 3D Plot
        figure(3)
        s = surf(X,Y,T);  hold on;
        s.FaceColor = 'interp';
        s = surf(X,-Y,T); hold off;
        s.FaceColor = 'interp';
        xlabel('x'); ylabel('y'); zlabel('T');
        set(gcf, 'Position',[1270,150,550,550]);
        view(45,24);            %view angle
        saveas(gcf,'Steady_3D.jpg')
        


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

        %% Plotting Step
        i=2;
        if (theta >= 0 && theta < 0.5) || (strcmp(TimeIntegrType, 'RungeKutta4'))
            i=50;  % raise the step because those 2 methods need too small Δt
        end

        %% Critical Δt
        % Plot
        if (theta >= 0 && theta < 0.5)
            figure(1)
            plot(dt_crit(1,:),dt_crit(2,:));   hold on;
            plot(dimX*dimY,dt,'o');            hold off;
            set(gcf, 'Position',[50,100,1500,850]);
            xlabel('Number of grid points'); ylabel('Δt [s]'); title('Critical Δt');
            text(dimX*dimY+25,dt+dt/13, 'Δt');

            % Display text
            format long
            if dt>dt_crit(2,i_min)
                fprintf(2, '\nWarning! Δt is too high!\n \tRecommended value, approx.:\t %s \n', num2str(dt_crit(2,i_min),3))
            else
                fprintf('Recommended Δt, approx.:\t %s \n', num2str(dt_crit(2,i_min),3))
            end
        end

        %% 3D Plot
        filename = 'Unsteady_3D.gif';
        figure(2)
        for t=1.0:i:size(T,3)
            s = surf(X,Y,T(:,:,t));        hold on;
            s.FaceColor = 'interp';   % values (colors) on nodes and not on elements
            s = surf(X,-Y,T(:,:,t));       hold off;
            s.FaceColor = 'interp';
            xlabel('x'); ylabel('y'); zlabel('T');
            title("t = " + t*dt + "[s]")
            set(gcf, 'Position',[10,150,550,550]);
            view(45,24);            %view angle

            % Record into a GIF
            drawnow
            frame= getframe(gcf);
            im= frame2im(frame);
            [imind,cm] = rgb2ind(im,64);
            if t == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
            end
        end

        %% Contour plot
        filename = 'Unsteady_Contour.gif';
        figure(4)
        for t=1.0:i:size(T,3)
            contour(X,Y,T(:,:,t));      hold on;
            contour(X,-Y,T(:,:,t));     hold on;
            plot(X(1,:),Y(1,:),'k');    hold on;
            plot(X(1,:),-Y(1,:),'k');   hold off;
            text(l-l/6,-Y(1,1),['t = ', num2str(t*dt)],'fontweight','bold')
            xlabel('x'); ylabel('y');
            if strcmp(shape, 'linear')
                axis([0 l -h1/2 h1/2 -1 1]) %limits of x & y axes
            else
                axis([0 l -h1/2*1.15 h1/2*1.15 -1 1]) %limits of x & y axes
            end
            hcb=colorbar;
            set(get(hcb,'Title'),'String','T');     % insert label above colorbar
            set(gcf, 'Position',[570,150,620,550]);

            %Record into a GIF
            drawnow
            frame= getframe(gcf);
            im= frame2im(frame);
            [imind,cm] = rgb2ind(im,64);
            if t == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
            end
        end

        %% 2D Plot
        filename = 'Unsteady_2D.gif';
        figure(3)
        for t=1.0:i:size(T,3)
            s = pcolor(X,Y,T(:,:,t));        hold on;
            s.FaceColor = 'interp';
            s = pcolor(X,-Y,T(:,:,t));       hold off;
            s.FaceColor = 'interp';
            text(l-l/6,-Y(1,1),['t = ', num2str(t*dt)],'fontweight','bold')
            xlabel('x'); ylabel('y');
            if strcmp(shape, 'linear')
                axis([0 l -h1/2 h1/2 -1 1]) %limits of x & y axes
            else
                axis([0 l -h1/2*1.15 h1/2*1.15 -1 1]) %limits of x & y axes
            end
            hcb=colorbar;
            set(get(hcb,'Title'),'String','T');
            set(gcf, 'Position',[1200,150,620,550]);

            %Record into a GIF
            drawnow
            frame= getframe(gcf);
            im= frame2im(frame);
            [imind,cm] = rgb2ind(im,64);
            if t == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
            end
        end
end