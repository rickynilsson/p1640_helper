function template_p1640w = IRTFtoP1640werr(irtf_template_file,wavelengths)

%IRTFTOP1640WERR: Get corresponding P1640 fluxes of IRTF spectral calibrator stars
% The wavelength solution varies slightly for diffent PCXP extractions. It is recommended that you find the most accurate wavelength solution before using this code to obtain your template spectrum for SRF calculation. If a wavelength scale is not supplied to function call, the default SPECTRUM.LAMBDA below will be used.  
% SPECTRUM.LAMBDA = linspace(995,1769,32) is the wavelength range as set by edges of filter transmission function (where transmission drops to 50%)
%
% Syntax: template_p1640w = IRTFtoP1640werr(irtf_template_file,wavelengths)
%
% Inputs:
%   irtf_template_file: Full path and filename of IRTF calibrator
%   wavelengths (optional): Wavelength scale as column vector
%
% Outputs:
%   template_p1640w: 32 x 3 array with [wavelength flux flux_err]
%
% Example: 
%   output_template_p1640w = IRTFtoP1640('/IRTF_091201/G5V_HD165185.fits',linspace(995,1769,32))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files or other data files required: Spectral template file from IRTF Spectral Library (http://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/)
%
% See also: N/A

% Author: Ricky Nilsson, Ph.D., Astronomy
% Department of Astronomy, California Institute of Technology, CA, USA
% email address: ricky@caltech.edu  
% Website: http://www.caltech.edu/~rnilsson
% November 2017; Last revision: 27-Nov-2017

%------------- BEGIN CODE --------------

SPECTRUM.SPECCAL = irtf_template_file;

if ~exist('wavelengths','var')
    SPECTRUM.LAMBDA = linspace(950,1750,32); % default wavel
    print('No wavelengths supplied, using default range and scale ...')
else
    SPECTRUM.LAMBDA = wavelengths; % supplied wavel
end

SPECTRUM.USE_BANDS = [1:32];

slashind = findstr(SPECTRUM.SPECCAL,'/');
TEMPLATENAME = SPECTRUM.SPECCAL(slashind(end)+1:end);

templatestar_spect = fits_read([SPECTRUM.SPECCAL]); % Read in calibrator template. IRTF data is in 3 rows with lambda, F, and dF
templatestar_spect = [1e3*templatestar_spect(1,:)' 1e3*templatestar_spect(2,:)' 1e3*templatestar_spect(3,:)']; % convert to column vectors, wavelength from micron to nm, and flux from W/m^2/nm to erg/s/cm^2/nm
templatestar_spect(any(isnan(templatestar_spect),2),:) = [];

w_instr = SPECTRUM.LAMBDA; % wavelength scale of instrument (nm)
dw_instr = diff(w_instr(1:2))*1.939; % The effective FWHM of laser PSF after PCXP in Zimmerman et al. (2011) is ~70nm, but was measured by E. Bacchus to 1.939 channels for Bumpy PCXP

c = 299792.458; % speed of light (km/s)
v_r = 0; % radial velocity of template star (km/s)
w_tstar = templatestar_spect(:,1) / (1 + v_r / c); % wavelength axis corrected for any radial velocity (nm)
f_tstar = templatestar_spect(:,2); % flux axis (erg/s/cm^2/nm)
df_tstar = templatestar_spect(:,3); % flux errors (erg/s/cm^2/nm)


%w_tstar_instrrange_ind = find(w_tstar > w_instr(1) & w_tstar < w_instr(end));
%w_tstar_instrrange = w_tstar(w_tstar > w_instr(1) & w_tstar < w_instr(end));
%f_tstar_instrrange = f_tstar(w_tstar_instrrange_ind);

%%% Calculate error bars in each new bin %%%
df_tstar_p1640 = zeros(32,1);
for band = 1:32
	w_tstar_instrrange_ind = find((w_tstar >= w_instr(band)-dw_instr) & (w_tstar <= w_instr(band)+dw_instr)); % index range in irtf template for each bin/band
	df_tstar_p1640(band) = sqrt(sumsqr(df_tstar(w_tstar_instrrange_ind)));
end 

w_tstar_reg = w_tstar(1):0.1:w_tstar(end);
f_tstar_reg = interp1(w_tstar,f_tstar,w_tstar_reg,'spline','extrap');
figure()
stairs(w_tstar,f_tstar); xlim([w_instr(1) w_instr(end)])
hold on
stairs(w_tstar_reg,f_tstar_reg,'g-')


dw_tstar = 0.1; % resolution of template spectrum (on regular grid) in P1640 range (nm) 
xg = -2*dw_instr : dw_tstar : 2*dw_instr; % x axis from -2 x (resolution of instrument) to +2 x (resolution of instrument)
b = 0;
fwhm = dw_instr; % full-width-half-max of gaussian, same as spectral resolution of instrument (nm)
sigma = fwhm / (2 * sqrt(2 * log(2))); % sigma of gaussian
gaussf = exp(-((xg - b).^2)/(2 * sigma^2)); % gaussian y values

figure()
plot(w_tstar)
ylim([w_instr(1) w_instr(end)])
%hold on
%plot(w_tstar_instrrange_reg,'g-')
title('Check continuity of IRTF template star wavelength scale in P1640 range. Will be interpolated onto regular grid.')
tstar_ldiff = max(diff(w_tstar(w_tstar > w_instr(1) & w_tstar < w_instr(end))));
if tstar_ldiff > max(diff(w_instr))
	disp('Spectrum has missing values at P1640 wavelengths! Do not use!')
	template_p1640w = [];
	return
end

figure()
stairs(w_tstar,f_tstar)
hold on

f_tstar_degraded = conv(f_tstar_reg,gaussf,'same');
%f_tstar_degraded(isnan(f_tstar_degraded)) = 0;
scalefactor_for_degraded = sum(f_tstar_reg)/sum(f_tstar_degraded);
stairs(w_tstar_reg,f_tstar_degraded * scalefactor_for_degraded,'g')

xlim([min(w_instr) max(w_instr)])

templatestar_spect_degraded = [w_tstar_reg' (f_tstar_degraded * scalefactor_for_degraded)'];
%f_tstar_degraded_p1640w = resample(templatestar_spect_degraded(:,2),round(dw_tstar*1000),(round(dw_instr*1000)));
f_tstar_degraded_p1640w = interp1(w_tstar_reg,templatestar_spect_degraded(:,2),w_instr,'linear','extrap');
f_tstar_p1640w = interp1(w_tstar,f_tstar,w_instr,'linear','extrap');

templatestar_spect_p1640w = [w_instr; f_tstar_degraded_p1640w];

plot(w_instr,f_tstar_degraded_p1640w,'k+')
plot(w_instr,f_tstar_p1640w,'r+')

ylabel('Flux (erg/s/cm^2/nm)')
xlabel('Wavelength (nm)')
title(TEMPLATENAME,'Interpreter','none')

legend('Full resolution template spectrum','Degraded template spectrum','P1640 flux from interpolation in degraded template spectrum','P1640 flux from interpolation in full-res template spectrum')

template_p1640w = [w_instr' f_tstar_degraded_p1640w' df_tstar_p1640]

end

%------------- END OF CODE --------------