if strcmpi(yeartype, 'Cal BP') == 1 || strcmpi(yeartype, 'CalBP') == 1
	yearlabel = 'cal yr BP';
elseif strcmpi(yeartype, 'BCE/CE') == 1
	yearlabel = 'cal yr BCE/CE';
else
	error('Please specify a valid year type (see help for options)')
end

% set plot settings depending on if there is a reservoir correction or not
if strcmp(curvetype, 'atm') == 1
	if isnan(resage) == 0 || isnan(reserr) == 0
		extralabel = 1; % extra line of text detailing reservoir correction info
		plot14Cextra = 1; % plot both original 14 determination & reservoir corrected one
		if strcmpi(dettype,'14cyr') == 1
			reslabel = 'R(t)';
		elseif strcmpi(dettype,'f14c') == 1
			reslabel = '\delta^1^4R';
		end
	else
		extralabel = 0;
		plot14Cextra = 0;
	end
elseif strcmp(curvetype, 'mar') == 1
	extralabel = 1;
	plot14Cextra = 1;
	if strcmpi(dettype,'14cyr') == 1
		reslabel = '\DeltaR';
	elseif strcmpi(dettype,'f14c') == 1
		reslabel = '\delta^1^4R';
	end
end


% prep age range for plot window
calprob2 = calprob(cumsum(calprob(:,2)) > 0.001 & cumsum(calprob(:,2)) < 0.999);
if strcmpi(yeartype,'BCE/CE') == 1
	calprob2 = flipud(calprob2);
	curvecal = (curvecal-1950) * -1;
end
yrrng = (calprob2(end,1) - calprob2(1,1))/2;
syr = (10^2) * round((calprob2(1,1)-yrrng) / (10^2)); % round to nearest hundred for nice plot limits
eyr = (10^2) * round((calprob2(end,1)+yrrng) / (10^2));
ind = find(curvecal >= syr & curvecal <= eyr);
curvecal = curvecal(ind);
curve14c = curve14c(ind);
curve14cerr = curve14cerr(ind);

if exist('OCTAVE_VERSION', 'builtin') ~= 0 % very simple output plot for octave users
	
	figure(14)
	clf
	plot(calprob(:,1),calprob(:,2),'b-')
	xlim([syr eyr])
	xlabel(['Age (',yearlabel,')'])
	ylabel('Calbriated age probability')
	title('Simplified output plot for OCTAVE')
	
else % otherwise continue with more fancy plot for matlab users
	
	figure(14)
	clf
	
	%----- Plot ProbDistFunc
	axpdf = axes;
	axes(axpdf)
	area(calprob(:,1),calprob(:,2),'edgecolor','none')
	axpdfylims = ylim;
	% axpdfxlims = xlim; % not used
	area(calprob(:,1),calprob(:,2)*0.2,'edgecolor',[0 0 0],'facecolor',[0.9 0.9 0.9])
	hold on
	for i = 1:size(p95_4,1)
		if strcmpi(yeartype,'Cal BP') == 1 || strcmpi(yeartype,'CalBP') == 1
			area( calprob(calprob(:,1) <= p95_4(i,1) & calprob(:,1) >= p95_4(i,2),1)  , calprob(calprob(:,1) <= p95_4(i,1) & calprob(:,1) >= p95_4(i,2),2)*0.2,'edgecolor','none','facecolor',[0.56 0.56 0.66])
		elseif strcmpi(yeartype,'BCE/CE') == 1
			area( calprob(calprob(:,1) >= p95_4(i,1) & calprob(:,1) <= p95_4(i,2),1)  , calprob(calprob(:,1) >= p95_4(i,1) & calprob(:,1) <= p95_4(i,2),2)*0.2,'edgecolor','none','facecolor',[0.56 0.56 0.66])
		end
	end
	for i = 1:size(p68_2,1)
		if strcmpi(yeartype,'Cal BP') == 1 || strcmpi(yeartype,'CalBP') == 1
			area( calprob(calprob(:,1) <= p68_2(i,1) & calprob(:,1) >= p68_2(i,2),1)  , calprob(calprob(:,1) <= p68_2(i,1) & calprob(:,1) >= p68_2(i,2),2)*0.2,'edgecolor','none','facecolor',[0.5 0.5 0.6])
		elseif strcmpi(yeartype,'BCE/CE') == 1
			area( calprob(calprob(:,1) >= p68_2(i,1) & calprob(:,1) <= p68_2(i,2),1)  , calprob(calprob(:,1) >= p68_2(i,1) & calprob(:,1) <= p68_2(i,2),2)*0.2,'edgecolor','none','facecolor',[0.5 0.5 0.6])
		end
	end
	
	%----- Plot cal curve
	axcurve = axes;
	axes(axcurve)
	xdata = curvecal;
	ydata = curve14c;
	onesig = curve14cerr;
	if strcmpi(dettype,'f14c') == 1
		ydata = exp(ydata/-8033);
		onesig = ydata.*onesig/8033;
	end
	fill([xdata' fliplr(xdata')],[ydata'+2*onesig' fliplr(ydata'-2*onesig')],[0.8 0.8 0.8],'edgecolor','none');
	hold on
	fill([xdata' fliplr(xdata')],[ydata'+onesig' fliplr(ydata'-onesig')],[0.6 0.6 0.6],'edgecolor','none');
	axcurveylims = ylim;
	axcurvexlims = [min(curvecal) max(curvecal)];
	
	
	%----- Plot raw data if intcal13 or intcal20 is selected
	if strcmpi('intcal13',calcurve) == 1 || strcmpi('intcal20',calcurve) == 1
		
		axraw = axes;
		axes(axraw)
		
		if strcmpi('intcal13',calcurve) == 1
			
			rd = load('IntCal13 raw data.txt');
			% now trim out the bits that weren't actually used in IntCal13
			rd_trees = rd(rd(:,1) >= 1 & rd(:,1) <= 8, :);
			rd_other = rd(rd(:,1) >= 9, :);
			rd_trees = rd_trees(rd_trees(:,3) <= 13900, :);
			rd_other = rd_other(rd_other(:,3) >= 13900, :);
			rd = [rd_trees; rd_other];
			if strcmpi(yeartype,'BCE/CE') == 1
				rd(:,3) = (rd(:,3)-1950) * -1;
				ind = find(rd(:,3) <= curvecal(1) & rd(:,3) >= curvecal(end));
			else
				ind = find(rd(:,3) >= curvecal(1) & rd(:,3) <= curvecal(end));
			end
			
			raw_cal = rd(ind,3);
			raw_calsigma = rd(ind,5);
			raw_14c = rd(ind,6);
			raw_14csigma = rd(ind,7);
			if strcmpi(dettype,'f14c') == 1
				raw_14c = exp(raw_14c/-8033);
				raw_14csigma = raw_14c.*raw_14csigma/8033;
			end
			
		elseif strcmpi('intcal20',calcurve) == 1
			
			rd = load('IntCal20 raw data.txt');
			if strcmpi(yeartype,'BCE/CE') == 1
				rd(:,3) = (rd(:,3)-1950) * -1;
				ind = find(rd(:,3) <= curvecal(1) & rd(:,3) >= curvecal(end));
			else
				ind = find(rd(:,3) >= curvecal(1) & rd(:,3) <= curvecal(end));
			end
			
			raw_cal = rd(ind,3);
			raw_calsigma = rd(ind,4);
			raw_14c = rd(ind,5);
			raw_14csigma = rd(ind,6);
			if strcmpi(dettype,'f14c') == 1
				raw_14c = exp(raw_14c/-8033);
				raw_14csigma = raw_14c.*raw_14csigma/8033;
			end
			
		end
		
		hbars = sbars(raw_cal, raw_14c, raw_calsigma, raw_14csigma);
		set(hbars,'color',[132 193 150]/256);
		
		axrawylims = [min(raw_14c-raw_14csigma) max(raw_14c+raw_14csigma)];
		% axrawxlims = xlim;  % not used
		
	end
	
	%----- Plot 14C age normal distribution(s)
	axgauss = axes;
	axes(axgauss);
	
	if strcmpi(dettype,'14cyr') == 1
		gaussrange = (labdet-4*laberr:labdet+4*laberr);
		if plot14Cextra == 1
			gaussrangeorig = (labdetorig-4*laberrorig:labdetorig+4*laberrorig);
		end
	elseif strcmpi(dettype,'f14c') == 1
		gaussrange = (labdet-4*laberr:0.00001:labdet+4*laberr);
		if plot14Cextra == 1
			gaussrangeorig = (labdetorig-4*laberrorig:0.00001:labdetorig+4*laberrorig);
		end
	end
	gauss = normpdf(gaussrange,labdet,laberr);
	
	patch(gauss, gaussrange, 'blue');
	% axgaussylims = ylim; % not used
	axgaussxlims = xlim;
	axgaussxlims(2) = axgaussxlims(2)*5;
	if plot14Cextra == 1
		gaussorig = normpdf(gaussrangeorig, labdetorig, laberrorig);
		a = patch(gaussorig, gaussrangeorig, 'blue');
		set(a,'edgecolor','none','facecolor',[0.8 0.5 0.5]);
		hold on
	end
	a = patch(gauss, gaussrange, 'blue');
	set(a,'edgecolor',[0 0 0],'facecolor',[0.7 0.4 0.4]);
	
	
	%----- set plot settings by axis, starting from back layer to front layer
	
	axes(axcurve)
	xlim(axcurvexlims)
	if revxdir == 1
		set(gca, 'XDir', 'reverse')
	end
	set(gca,'color','none')
	if strcmpi(dettype,'14cyr') == 1
		lab1 = ylabel('Conventional ^1^4C age (^1^4C yr BP)');
	elseif strcmpi(dettype,'f14c') == 1
		lab1 = ylabel('^1^4C activity (F^1^4C)');
	end
	lab2 = xlabel(['Calibrated age (',yearlabel,')']);
	set( gca, 'TickDir', 'out' );
	if strcmpi('intcal13',calcurve) == 1 || strcmpi('intcal20',calcurve) == 1
		ylim(axrawylims)
	else
		ylim(axcurveylims)
	end
	yt=get(gca,'ytick');
	if strcmpi(dettype,'14cyr') == 1
		ytl=textscan(sprintf('%1.0f \n',yt),'%s','delimiter','');
		set(gca,'yticklabel',ytl{1})
	end
	xt=get(gca,'xtick');
	xtl=textscan(sprintf('%1.0f \n',xt),'%s','delimiter','');
	set(gca,'xticklabel',xtl{1})
	
	axes(axpdf)
	set(gca,'color','none')
	xlim(axcurvexlims)
	ylim(axpdfylims)
	set(gca,'xticklabel',[]);
	set(gca,'xtick',[]);
	set(gca,'yticklabel',[]);
	if revxdir == 1
		set(gca, 'XDir', 'reverse')
	end
	set(gca,'ytick',[]);
	
	axes(axgauss)
	set(gca,'color','none')
	xlim(axgaussxlims)
	set(gca,'xticklabel',[]);
	set(gca,'xtick',[]);
	set(gca,'yticklabel',[]);
	set(gca,'ytick',[]);
	if strcmpi('intcal13',calcurve) == 1 || strcmpi('intcal20',calcurve) == 1
		ylim(axrawylims)
	else
		ylim(axcurveylims)
	end
	
	axes(axcurve) % bring curve forward
	
	if strcmpi('intcal13',calcurve) == 1 || strcmpi('intcal20',calcurve) == 1
		axes(axraw)
		hold on
		set(gca,'color','none')
		xlim(axcurvexlims)
		ylim(axrawylims)
		set(gca,'xticklabel',[]);
		set(gca,'xtick',[]);
		set(gca,'yticklabel',[]);
		set(gca,'ytick',[]);
		if revxdir == 1
			set(gca, 'XDir', 'reverse')
		end
	else
		axes(axpdf)
	end
	
	
	
	%----- Plot some text on the final axis
	
	% Top left text box: MatCal version and cal curve used
	verboxstr = [matcalvers,newline,calcurve, ' ', cite];
	vbh = annotation('textbox',get(gca,'position'),'String',verboxstr);
	set(vbh, 'linestyle','none')
	set(vbh, 'horizontalalignment','left')
	set(vbh, 'verticalalignment','top')
	
	% Top right text box: Raw and calibrated age info
	if strcmpi(dettype,'14cyr') == 1
		ageboxstr = ['^1^4C det.: ',num2str(labdetorig),' \pm ',num2str(laberrorig),' ^1^4C yr BP'];
	elseif strcmpi(dettype,'f14c') == 1
		ageboxstr = ['^1^4C det.: ',num2str(labdetorig,'%.5f'),' \pm ',num2str(laberrorig,'%.5f'),' F^1^4C'];
	end
	
	if extralabel == 1
		if strcmpi(dettype,'14cyr') == 1
			ageboxstr = [ageboxstr,newline,reslabel,': ',num2str(resage),' \pm ',num2str(reserr), ' ^1^4C yr'];
		elseif strcmpi(dettype,'f14c') == 1
			ageboxstr = [ageboxstr,newline,reslabel,': ',num2str(resage),' \pm ',num2str(reserr)];
		end
	end
	
	ageboxstr = [ageboxstr,newline];
	
	if size(p95_4,1) == 1
		ageboxstr = [ageboxstr,newline,'Cal age 95.45% HPD interval:'];
	else
		ageboxstr = [ageboxstr,newline,'Cal age 95.45% HPD intervals:'];
	end
	for i = 1:size(p95_4,1)
		ageboxstr = [ageboxstr,newline,num2str(floor(p95_4(i,3)*1000)/10),'%: ',num2str(p95_4(i,1)),' to ',num2str(p95_4(i,2)),' ',yearlabel];
	end
	
	abh = annotation('textbox',get(gca,'position'),'String',ageboxstr);
	set(abh, 'linestyle','none')
	set(abh, 'horizontalalignment','right')
	set(abh, 'verticalalignment','top')
	
	
	%----- Warning if tail of 14C date or cal age is near limit of cal curve
	if warn == 1
		title('Warning! Age calibration may exceed limits of calibration curve.')
	end
	
	%----- Uniform fonts and appearance
	
	set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
	set(lab1,'FontWeight','bold')
	set(lab2,'FontWeight','bold')
	set(gcf,'color',[1 1 1]);
	
end

%----- Prep plot for export
if printme == 1
	
	% set figure size (cm)
	xSize = plotsize;
	ySize = plotsize;
	
	% set paper size (cm)f
	set(gcf,'PaperUnits','centimeters')
	Y = plotsize+2;
	X = plotsize+2;
	set(gcf, 'PaperSize',[X Y])
	
	% put figure in centre of paper
	xLeft = (X-xSize)/2;
	yBottom = (Y-ySize)/2;
	set(gcf,'PaperPosition',[xLeft yBottom xSize ySize])
	% make background white
	set(gcf,'InvertHardcopy','on');
	set(gcf,'color',[1 1 1]);
	
	print(figure(14), '-dpdf', '-painters', ['MatCal ',num2str(labdetorig),char(177),num2str(laberrorig),'.pdf']);
end
