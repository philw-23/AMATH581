function plotSolutions(t_span, sol, x_grid, y_grid, z_grid, N, p, filename)
    
    % Plot GIFs
    if strcmp(p, 'gif')
        filename = strcat(filename, '.gif'); % add .gif to file name
        for i=1:length(t_span) % iterate through timespan
            y_plot = reshape(sol(i, :), N, N, []); % get value for time
            y_plot = real(ifftn(y_plot)); % get real value in real space
            isosurface(x_grid, y_grid, z_grid, y_plot); % plot surface
            xlabel('x'); ylabel('y'); zlabel('z'); % Label axes
            title(['Phi at time t = ', num2str(t_span(i)), ' seconds'])
            
            % Get frame and convert to image
            frame = getframe(gcf); % get the frame for f
            image = frame2im(frame); % convert frame to image
            [imind, cmap] = rgb2ind(image, 256); % convert to indexed image
    
            % Write out to gif file
            if i == 1 % generate initial image write
                imwrite(imind, cmap, filename, 'gif', 'LoopCount', inf);
            else % add additional frames
                imwrite(imind, cmap, filename, 'gif', 'WriteMode', 'append');
            end  
        end
        clf; % clear figure
    end

    % Plot time progressions (not working)
    if strcmp(p, 'prog')
        filename = strcat(filename, ' Progression'); % Add piece to title
        f = figure(); % default figure
        for i=1:length(t_span)
            subplot(3, 3, i)
            y_plot = reshape(sol(i, :), N, N, []); % get value for time
            y_plot = real(ifftn(y_plot)); % get real value in real space
            isosurface(x_grid, y_grid, z_grid, y_plot); % plot surface
            title(['t = ', num2str(t_span(i)), ' seconds']); % title
            view(3); % set view dimension
            camlight; lighting gouraud; % set lighting defaults
        end
        saveas(f, strcat(filename, '.jpg')); % save as jpeg
        clf; % clear figure
    end

end