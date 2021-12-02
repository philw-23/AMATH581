% Plot ode45 solution parameters
function plotSolutions(t_span, sol, N, x_grid, y_grid, method)
    
    for j=1:length(t_span)
        
        uplot = reshape(sol(j,1:N^2),N,[]); % u values
        vplot = reshape(sol(j,N^2+1:end),N,[]); % v values
        
        % convert out of fourier space if solution is from FFT
        if strcmp(method, 'fft') % FFT solutions
            uplot = real(ifft2(uplot));
            vplot = real(ifft2(vplot));
        end

        % Plot u
        subplot(1, 2, 1); % Left Subplot
        pcolor(x_grid,y_grid,uplot);
        shading interp; axis square;
        title('u');

        % Plot v
        subplot(1, 2, 2); % Right Subplot
        pcolor(x_grid,y_grid,vplot);
        shading interp; axis square; % no grid lines and square axes
        title('v');
        
        % Draw plot
        drawnow;
        pause(0.1);

    end
end