function  [p95_4, p68_2, calprob, medage] = matcal(labdet, laberr, calcurve, yeartype, varargin)
% [p95_4, p68_2, calprob, medage] = matcal(labdet, laberr, calcurve, yeartype)
%
% Function for 14C age calibration using Bayesian higher posterior
% density analysis of a probability density function of calibrated age.
%
% Please see manuscript for more detailed information:
% Lougheed, B.C. & Obrochta, S.P. (2016). MatCal: Open Source Bayesian
% 14C Age Calibration in Matlab. Journal of Open Research Software. 4(1),
% p.e42. DOI: http://doi.org/10.5334/jors.130
%
% Required input parameters
% =================================
%
% labdet:    Lab 14C determination (default unit is 14C yr BP).
%
% laberr:    Lab 14C determination uncertainty (1 sigma).
%
% calcurve:  String specifying calibration curve to use, select from
%            the following (not case sensitive):
%            'IntCal20', 'Marine20', 'SHCal20', 'IntCal13', 'Marine13',
%			 'SHCal13, 'IntCal09', 'Marine09', 'IntCal04', 'Marine04',
%			 'SHCal04, 'IntCal98', 'Marine98'
%
% yeartype:  String specifying how to report calibrated age.
%            Choices are 'CalBP' or 'BCE/CE'. (Not case sensitive)
%
% Optional input parameters
% =================================
%
% resage:    Optional (parameter name and value). Specify reservoir
%            age in 14C yr. Reservoir age is R(t) in the case of an 
%			 atmospheric calibration curve, and DeltaR in the case of
%			 a marine curve. (default = 0)
%            e.g. 'resage',320 for a reservoir age of 320
%
% reserr:    Optional (parameter name and value). Specify a 1 sigma
%            uncertainty for your chosen resage (default = 0)
%            e.g. 'reserr',50 for an uncertainty of 50
%
% dettype:   Optional (parameter name and value). Choose to handle
%            lab determination and plot as F14C or 14C years (default).
%            Options are '14cyr' (default) or 'f14c'. Not case sensitive.
%            e.g. 'dettype','f14c' to activate F14C mode.
%            If F14C is chosen then lab determination (labdet, laberr)
%            must be inputted as F14C, with resage (and reserr)
%            inputted as the F14C depletion ratio relative to the chosen
%            calibration curve, in per mil notation. This notation has been
%            referred to as d14R by Soulet et al (2016), doi: 10.1017/RDC.2015.22.            
%
% plot:      Optional (parameter name and value). Return a calibration
%            plot to Figure 14. The plot displays the 1 and 2
%            sigma ranges of the calibration curve. The calibraiton
%            curve raw data is also shown if IntCal13 is selected.
%            Specify 1 to plot and 0 not to plot. (default = 1 for Matlab
%            users; default = 0 for Octave users)
%            e.g. 'plot',0 not to plot
%
% saveplot:  Optional (parameter name and value). Save an Adobe PDF of
%            the calibration plot to your working directory. Specify 1
%            to save and 0 not to save. (default = 0) Will be ignored
%            if plotting has been disabled.
%            e.g. 'saveplot',1 to save to your working directory.
%
% plotsize:  Optional (parameter name and value). Set the width and height
%            of the printed figure in cm. (default = 16).
%            e.g. 'plotsize',10 for 10 cm.
%
% fontsize:  Optional (parameter name and value). Set the value of the font
%            size in the output plot. (default = 8)
%            e.g. 'fontsize',12 for a font size of 12.
%
% revxdir:   Optional (parameter name and value). Reverse the plot x-axis.
%            Specify 1 to reverse and 0 not to reverse. (default = 0)
%            e.g. 'revxdir',1 to reverse the x-axis.
%
% Output data
% =================================
%
% p95_4:     n by 3 matrix containing 95.45% calibrated age probability
%            range interval(s) calculated using highest posterior density.
%            Each row contains a probability range in Cols 1 and 2, and
%            the associated probability for that range in Col 3.
%            Probabilities are normalised to between zero and one.
%
% p68_2:     Same as p95_4, but for the 68.27% calibrated range.
%
% calprob:   n by 2 matrix containing an annualised calibrated age
%            probability density function for implementation in, e.g.,
%            age modelling. n is the annualised length of the chosen
%            calibration curve. Col 1 is a series of annual cal ages,
%            Col 2 contains their associated probability. All probabilities
%            are normalised such that they sum to 1.
%
% medage:    Median age calculated from calprob.
%
% Functional examples
% =================================
%
% [p95_4, p68_2, prob, medage] = matcal(1175, 30, 'IntCal20', 'BCE/CE');
% Calibrate a 14C age of 1175±30 14C yr BP using IntCal20 with output in BCE/CE.
%
% [p95_4, p68_2, prob, medage] = matcal(23175, 60, 'Marine20', 'CalBP',...
% 'resage', -50, 'reserr', 100, 'saveplot', 1);
% Calibrate a 14C age of 23175±60 14C yr BP using Marine20, with output in
% Cal BP, with delta-R of -50±100 14C yr and save a copy of the plot to
% your working directory as an Adobe PDF.
%
% [p95_4, p68_2, prob, medage] = matcal(1175, 30, 'IntCal20', 'CalBP', 'plot', 0);
% Calibrate a 14C age of 1175±50 14C yr BP using IntCal20, with output in
% Cal BP and disable the plot window.
%
% ------------
%
% MatCal 3.1 (2021-05-18)
% Originally written using Matlab 2012a, tested compatible with 2019a.
% No toolboxes required.
% Please see manuscript for more information:
% http://doi.org/10.5334/jors.130
matcalvers = 'MatCal (Lougheed and Obrochta, 2016), ver 3.1';

if nargin < 4
	error('Not enough input parameters (see help for instructions)')
end

% Optional parameters input parser (matlab's terrible parse of varargin)
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.FunctionName='matcal';
defaultresage = NaN;
defaultreserr = NaN;
defaultplotme = 1;
defaultprintme = 0;
defaultplotsize = 16;
defaultfontsize = 8;
defaultrevxdir = 0;
defaultdettype = '14cyr';
if exist('OCTAVE_VERSION', 'builtin') ~= 0
	addParamValue(p,'resage',defaultresage,@isnumeric); %#ok<*NVREPL>
	addParamValue(p,'reserr',defaultreserr,@isnumeric);
	addParamValue(p,'plot',defaultplotme,@isnumeric); 
	addParamValue(p,'saveplot',defaultprintme,@isnumeric);
	addParamValue(p,'plotsize',defaultplotsize,@isnumeric);
	addParamValue(p,'fontsize',defaultfontsize,@isnumeric);
	addParamValue(p,'revxdir',defaultrevxdir,@isnumeric);
	addParamValue(p,'dettype',defaultdettype);
else
	if datenum(version('-date'))>datenum('May 19, 2013')
		addParameter(p,'resage',defaultresage,@isnumeric);
		addParameter(p,'reserr',defaultreserr,@isnumeric);
		addParameter(p,'plot',defaultplotme,@isnumeric);
		addParameter(p,'saveplot',defaultprintme,@isnumeric);
		addParameter(p,'plotsize',defaultplotsize,@isnumeric);
		addParameter(p,'fontsize',defaultfontsize,@isnumeric);
		addParameter(p,'revxdir',defaultrevxdir,@isnumeric);
		addParameter(p,'dettype',defaultdettype);
	else
		addParamValue(p,'resage',defaultresage,@isnumeric);
		addParamValue(p,'reserr',defaultreserr,@isnumeric);
		addParamValue(p,'plot',defaultplotme,@isnumeric);
		addParamValue(p,'saveplot',defaultprintme,@isnumeric);
		addParamValue(p,'plotsize',defaultplotsize,@isnumeric);
		addParamValue(p,'fontsize',defaultfontsize,@isnumeric);
		addParamValue(p,'revxdir',defaultrevxdir,@isnumeric);
		addParamValue(p,'dettype',defaultdettype);
	end
end
parse(p,varargin{:});
resage = p.Results.resage;
reserr = p.Results.reserr;
plotme = p.Results.plot;
printme = p.Results.saveplot;
plotsize = p.Results.plotsize;
fontsize = p.Results.fontsize;
revxdir = p.Results.revxdir;
dettype = p.Results.dettype;

% Cal curve case and symbols
headerlines = 11;
if strcmpi(calcurve, 'IntCal20') == 1
	calcurve = 'IntCal20';
	cite = '(Reimer et al., 2020)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine20') == 1
	calcurve = 'Marine20';
	cite = '(Heaton et al., 2020)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'SHCal20') == 1
	calcurve = 'SHCal20';
	cite = '(Hogg et al., 2020)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'IntCal13') == 1
	calcurve = 'IntCal13';
	cite = '(Reimer et al., 2013)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine13') == 1
	calcurve = 'Marine13';
	cite = '(Reimer et al., 2013)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'SHCal13') == 1
	calcurve = 'SHCal13';
	cite = '(Hogg et al., 2013)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'IntCal09') == 1
	calcurve = 'IntCal09';
	cite = '(Reimer et al., 2009)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine09') == 1
	calcurve = 'Marine09';
	cite = '(Reimer et al., 2009)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'IntCal04') == 1
	calcurve = 'IntCal04';
	cite = '(Reimer et al., 2004)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine04') == 1
	calcurve = 'Marine04';
	cite = '(Hughen et al., 2004)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'SHCal04') == 1
	calcurve = 'SHCal04';
	cite = '(McCormac et al., 2004)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'IntCal98') == 1
	headerlines = 18;
	calcurve = 'IntCal98';
	cite = '(Stuiver et al., 1998)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine98') == 1
	headerlines = 18;
	calcurve = 'Marine98';
	cite = '(Stuiver et al., 1998)';
	curvetype = 'mar';
else
	error(['Calibration curve "',calcurve,'" unknown. Please specify a valid calibration curve (see help for options)'])
end

% store original 14C ages in memory and process reservoir age
if isnan(resage) == 1
	resage = 0;
end
if isnan(reserr) == 1 
	reserr = 0;
end
if strcmpi(dettype,'14cyr') == 1
	labdetorig = labdet;
	laberrorig = laberr;
	labdet = labdet - resage;
	laberr = sqrt(laberr^2 + reserr^2);
	f14cdet = exp(labdet/-8033);
	f14cerr = f14cdet*laberr/8033;
elseif strcmpi(dettype,'f14c') == 1
	labdetorig = labdet;
	laberrorig = laberr;
	labdet = labdet / (resage/1000+1); % d14R to F14R
	laberr = sqrt(laberr^2 + (reserr/1000)^2); % d14R to F14R
	f14cdet = labdet;
	f14cerr = laberr;
end

if strcmpi(yeartype, 'Cal BP') ~= 1 && strcmpi(yeartype, 'CalBP') ~= 1 && strcmpi(yeartype, 'BCE/CE') ~= 1
	error('Please specify a valid year type (see help for options)')
end

% open cal curve data
File = fopen([calcurve,'.14c']);
Contents = textscan(File,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(File);
curvecal = flipud(Contents{1});
curve14c = flipud(Contents{2});
curve14cerr = flipud(Contents{3});
curvef14 = exp(curve14c/-8033); % and also convert to F14 space
curvef14err = curvef14.*curve14cerr/8033; % and also convert to F14 space

% interpolate F14C cal curve to annual resolution
interpres = 1;
calprob(:,1) = curvecal(1):interpres:curvecal(end);
hicurvef14 = interp1(curvecal, curvef14, calprob(:,1));
hicurvef14err = interp1(curvecal, curvef14err, calprob(:,1));
% Calculate probability for every cal year in F14C space
% equation from e.g. p.261 in Bronk Ramsey, 2008. doi:10.1111/j.1475-4754.2008.00394.x
calprob(:,2) = exp(-((f14cdet - hicurvef14).^2)./(2 .* (f14cerr^2 + hicurvef14err.^2))) ./ ((f14cerr^2 + hicurvef14err.^2).^0.5) ;
calprob(:,2) = calprob(:,2) / sum(calprob(:,2)); % normalise to 1

warn = 0;
if strcmpi(dettype,'14cyr') == 1
	% throw warning if 4sigma of 14C age exceeds 14C age limits in cal curve
	if (labdet + 4*laberr) > max(curve14c) || (labdet - 4*laberr) < min(curve14c)
		warn = 1;
		warning(['4sigma range of 14C age ',num2str(labdetorig),char(177),num2str(laberrorig),' may exceed limits of calibration curve'])
	end
	% also throw warning if cal age PDF does not tail to zero at ends of cal curve (exceeds cal curve)
	if calprob(1,2) > 0.000001 || calprob(end,2) > 0.000001
		warn = 1;
		warning(['Calibrated age PDF for 14C age ',num2str(labdetorig),char(177),num2str(laberrorig),' may exceed limits of calibration curve'])
	end
elseif strcmpi(dettype,'f14c') == 1
	% throw warning if 4sigma of 14C age exceeds 14C age limits in cal curve
	if (labdet + 4*laberr) > max(hicurvef14) || (labdet - 4*laberr) < min(hicurvef14)
		warn = 1;
		warning(['4sigma range of F14C activity ',num2str(labdetorig),char(177),num2str(laberrorig),' may exceed limits of calibration curve'])
	end
	% also throw warning if cal age PDF does not tail to zero at ends of cal curve (exceeds cal curve)
	if calprob(1,2) > 0.000001 || calprob(end,2) > 0.000001
		warn = 1;
		warning(['Calibrated age PDF for F14C activity ',num2str(labdetorig),char(177),num2str(laberrorig),' may exceed limits of calibration curve'])
	end
end

% find 68.2% and 95.4% credible intervals using Bayesian highest posterior density (HPD)
% done by hpdcalc function in private folder
p68_2 = hpdcalc(calprob(:,1), calprob(:,2), 1, erf(1/sqrt(2)) );
p95_4 = hpdcalc(calprob(:,1), calprob(:,2), 1, erf(2/sqrt(2)) );

% calculate median (can't use interp1 because of potential repeat values)
[~, median_ind] = min(abs(cumsum(calprob(:,2))-0.5));
medage = round(mean(calprob(median_ind,1)));

% Convert output to BCE/CE if necessary
if strcmpi(yeartype,'BCE/CE') == 1
	medage = (medage-1950) * -1;
	calprob(:,1) = (calprob(:,1)-1950) * -1;
	p95_4(:,1:2) = (p95_4(:,1:2)-1950) * -1;
	p68_2(:,1:2) = (p68_2(:,1:2)-1950) * -1;
end

if plotme == 1
	% call plotting script from private folder
	matcalplot
end

end % end function
