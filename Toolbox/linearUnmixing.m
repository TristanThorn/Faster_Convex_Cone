function SO2 = linearUnmixing(PA, spectra_HbO2, spectra_Hb)
    % Ensure inputs are double precision for numeric stability
    PA = double(PA);
    spectra_HbO2 = double(spectra_HbO2);
    spectra_Hb = double(spectra_Hb);

    % Solve the constrained least squares problem enforcing non-negativity
    C = lsqnonneg([spectra_HbO2', spectra_Hb'], PA');

    if C(1) < 0
        SO2 = 0;
    elseif C(2) < 0
        SO2 = 1;
    else
        SO2 = C(1) / sum(C);
    end

    wavelengths = linspace(700, 900, 21);
    fitted_spectrum = C(1) * spectra_HbO2 + C(2) * spectra_Hb;

    figure;
    plot(PA, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Actual Spectrum');
    hold on;
    plot(fitted_spectrum, 'r*-', 'LineWidth', 1.5, 'DisplayName', 'Fitted Spectrum');

    xlabel('Wavelength (nm)');
    ylabel('PA Signal Intensity');
    title(sprintf('Linear Unmixing Fit (Predicted %.2f)', SO2));
    legend;
    set(gcf, 'Position', [200, 400, 600, 400]);
    grid on;
end
