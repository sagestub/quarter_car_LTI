function [h,x] = func_roadElevationProfile(k,L,dx,options)
    % [h,x] = func_roadElevationProfile(k,L,dx,options)
    %
    % Example : [h,x] = func_roadElevationProfile;
    %
    % reference: https://in.mathworks.com/matlabcentral/answers/511993-i-need-to-generate-standarad-iso-road-a-h-roughness-classification-from-3-to-9
    arguments
        k = [3:9]; % Constant 'k' Value as per Road roughness ISO 8608
        L = 2000; % Length Of Road Profile (m)
        dx = 10; % Sampling Interval (m)
        options.figure = true;
        options.fignum = {' '};
    end
    % Constant 'k' Value as per Road roughness ISO 8608
    k = reshape(k,1,[]);
    B = dx ; % Sampling Interval (m)
    dn = 1/L; % Frequency Band
    n0 = 0.1; % Spatial Frequency (cycles/m)
    N = round(L/dx); % Number of data points
    n = dn:dn:dn*N; % Spatial Frequency Band
    phi = 2*pi*rand(size(n)); % Random Phase Angle
    Amp = sqrt(dn)*(2.^k)'*(1e-3)*(n0./n);
    x = 0:B:L-B; % Abscissa Variable from 0 to L
    h = Amp * cos(2*pi*x'*n + phi)';
    if options.figure
        figure(options.fignum)
        hold on
        plot(x,h)
        xlabel('Distance (m)');
        ylabel('displacement Of PSD (m)');
        grid on
    end
end