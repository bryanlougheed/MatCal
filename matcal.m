function  [p95_4 p68_2 prob] = matcal(c14age, c14err, calcurve, yeartype, varargin)    
% [p68_2 p95_4 prob] = matcal(c14age, c14err, calcurve, yeartype, resage, reserr)
%
% Function for 14C age calibration using Bayesian statistical analysis of a
% probability density function of calibrated age.
%
% --- Input parameters ---
%
% c14age   = Lab 14C determination in 14C yr BP.
%
% c14err   = Lab 14C determination uncertainty (1 sigma) in 14C yr.
%
% calcurve = String specifying calibration curve to use, slect from:
%            'IntCal13', 'Marine13', 'SHCal13, 'IntCal09', 'Marine09'
%            'IntCal04', 'Marine04', 'SHCal04, 'IntCal98', 'Marine98'
%
% yeartype = String specifying how MatCal will report calibrated age.
%            Choices are 'Cal BP' or 'BCE/CE'.
%
% resage   = Optional. Specify reservoir age in 14C yr (default = 0).
%            R(t) in the case of atmospheric curve, delta-R in the case of
%            marine curve.
%            
% reserr   = Optional. Specify 1 sigma uncertainty for resage (default = 0)
%
% plot     = Optional. Default is to plot. Specify 0 for no plot
% 
% --- Output parameters ---
%
% p95_4    = n by 3 matrix containing 95.45% calibrated age probabilities
%            calculated using highest posterior density interval. Each row
%            contains a probability range in Cols 1 and 2, and the
%            associated probability for that range in Col 3. Probabilities
%            are normalised to between zero and one.
%
% p68_2    = Same as p68_2 but for the 68.27% calibrated age probabilities.
%
% prob     = n by 2 matrix containing a calibrated age probability
%            density function for implementation in, e.g., age modelling.
%            Col 1 is aseries of annual cal ages, Col 2 contains their
%            associated probability. Probabilities are normalised to
%            between zero and one.
%
% --- Plot ---
%
% A publication ready plot of the 14C age distribution, 14C age normal
% distribution with reservoir correction (if applicable), calibration curve
% (one sigma range) and probability density function of the calibrated age
% is displayed as figure(14).If IntCal13 is selected as calibration curve,
% the IntCal13 raw data is also shown in the background for reference. To
% save an Adobe PDF of the figure to your working directory, execute the
% following line in Matlab after running MatCal:
%
% print(figure(14), '-dpdf', 'name_of_your_figure.pdf');
%
% ------------
%
% Bryan C. Lougheed. Version 1.0 (May 2016). Written using MatLab 2012a.
% Feel free to modify for own use. Use of this script is at your own risk.

if nargin < 4
    error('Not enough input parameters (see help for instructions)')
end

%parse varargin

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.FunctionName='matcal';

defaultresage=0;
defaultreserr=0;
defaultprintme=1;

if datenum(version('-date'))<datenum('May 19, 2013')
    addParameter(p,'resage',defaultresage,@isnumeric);
    addParameter(p,'reserr',defaultreserr,@isnumeric);
    addParameter(p,'plot',defaultprintme,@isnumeric);
else
    addParamValue(p,'resage',defaultresage,@isnumeric);
    addParamValue(p,'reserr',defaultreserr,@isnumeric);
    addParamValue(p,'plot',defaultprintme,@isnumeric);
end

parse(p,varargin{:});
resage=p.Results.resage
reserr=p.Results.reserr
plotme=p.Results.plot

if resage==0
    extralabel=0;
    if reserr~=0
        error('Please specify reservoir age when specifying reservoir error')
    end
else
    extralabel=1;
end

% Cal curve case and symbols

headerlines = 11;

if strcmpi(calcurve, 'IntCal13') == 1
    calcurve = 'IntCal13';
    reslabel = 'R(t)';
    cite = '(Reimer et al., 2013)';
elseif strcmpi(calcurve, 'Marine13') == 1
    calcurve = 'Marine13';
    reslabel = '\DeltaR';
    cite = '(Reimer et al., 2013)';
    extralabel = 1;
elseif strcmpi(calcurve, 'SHCal13') == 1
    calcurve = 'SHCal13';
    reslabel = 'R(t)';
    cite = '(Hogg et al., 2013)';
elseif strcmpi(calcurve, 'IntCal09') == 1
    calcurve = 'IntCal09';
    reslabel = 'R(t)';
    cite = '(Reimer et al., 2009)';
elseif strcmpi(calcurve, 'Marine09') == 1
    calcurve = 'Marine09';
    reslabel = '\DeltaR';
    cite = '(Reimer et al., 2009)';
    extralabel = 1;
elseif strcmpi(calcurve, 'IntCal04') == 1
    calcurve = 'IntCal04';
    reslabel = 'R(t)';
    cite = '(Reimer et al., 2004)';
elseif strcmpi(calcurve, 'Marine04') == 1
    calcurve = 'Marine04';
    reslabel = '\DeltaR';
    cite = '(Hughen et al., 2004)';
    extralabel = 1;
elseif strcmpi(calcurve, 'SHCal04') == 1
    calcurve = 'SHCal04';
    reslabel = 'R(t)';
    cite = '(McCormac et al., 2004)';
elseif strcmpi(calcurve, 'IntCal98') == 1
    headerlines = 18;
    calcurve = 'IntCal98';
    reslabel = 'R(t)';
    cite = '(Stuiver et al., 1998)';
elseif strcmpi(calcurve, 'Marine98') == 1
    headerlines = 18;
    calcurve = 'Marine98';
    reslabel = '\DeltaR';
    cite = '(Stuiver et al., 1998)';
    extralabel = 1;
else
    error(['Please specify a valid calibration curve (see help for options)'])
end

if strcmpi(yeartype, 'Cal BP')
    yearlabel = 'Cal yr BP';
elseif strcmpi(yeartype, 'BCE/CE')
    yearlabel = 'BCE/CE';
else
    error(['Please specify a valid year type (see help for options)'])
end

% correct for reservoir age and store original ages in memory

c14ageorig = c14age;
c14errorig = c14err;

c14age = c14age - resage;
c14err = sqrt(c14err^2 + reserr^2);

% open cal curve data

File = fopen(['private/',calcurve,'.14c']);
Contents = textscan(File,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(File);
curvecal = flipud(Contents{1});
curve14c = flipud(Contents{2});
curve1sig = flipud(Contents{3});

% interpolate calibration curve to annual resolution

hicurvecal = curvecal(1):1:curvecal(end);
hicurve14c = interp1(curvecal, curve14c, hicurvecal);
hicurve1sig = interp1(curvecal, curve1sig, hicurvecal);

% Calculate probability for every cal year

prob = NaN(length(hicurvecal),2);

z = 0;
for i = 1:length(hicurvecal);

    z = z + 1;
    
    % equation from e.g. p.261 in Bronk Ramsey, 2008. doi:10.1111/j.1475-4754.2008.00394.x        
    % split equation into parts for sanity's sake 
    a = ( c14age - hicurve14c(i) )^2;
    b = 2 * (c14err^2 + hicurve1sig(i)^2);
    c = sqrt(c14err^2 + hicurve1sig(i)^2);
    prob(z,2) = exp(-a/b) / c;
    
    prob(z,1) = hicurvecal(i);
        
end

% normalise to 1

prob(:,2) = prob(:,2) / sum(prob(:,2)); 

% throw warning if PDF does not tail to zero on both sides (exceeds cal curve)
if prob(1,2) > 0.000001 || prob(end,2) > 0.000001
    warning(['Calibrated age for 14C age ',num2str(c14ageorig),'±',num2str(c14errorig),' may exceed limits of calibration curve'])
end

% find 1sig and 2sig intervals using highest posterior density (HPD)

hpd = prob(:,1:2);

hpd = sortrows(hpd, 2);

hpd(:,3) = cumsum(hpd(:,2));

% 1 sig

hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)), :);

hpd68_2 = sortrows(hpd68_2,1);

ind1 = find(diff(hpd68_2(:,1)) > 1);

if isempty(ind1) == 1;
    p68_2(1,1) = hpd68_2(end,1);
    p68_2(1,2) = hpd68_2(1,1);
    p68_2(1,3) = sum(hpd68_2(1:end,2));
else
    z = 0;
    for i = 1:length(ind1)
        z = z + 1;
        indy1(z,1) = ind1(i);
        z = z + 1;
        indy1(z,1) = ind1(i)+1;
    end 
    indy1 = [ 1 ; indy1; length(hpd68_2(:,1)) ];
    z=0;
    for i = 2:2:length(indy1)
        z = z+1;
        p68_2(z,1) = hpd68_2(indy1(i),1);
        p68_2(z,2) = hpd68_2(indy1(i-1),1);
        p68_2(z,3) = sum(hpd68_2(indy1(i-1):indy1(i),2));    
    end
    p68_2 = flipud(p68_2);
end


% 2 sig

hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)), :);

hpd95_4 = sortrows(hpd95_4,1);

ind2 = find(diff(hpd95_4(:,1)) > 1);

if isempty(ind2) == 1;
    p95_4(1,1) = hpd95_4(end,1);
    p95_4(1,2) = hpd95_4(1,1);
    p95_4(1,3) = sum(hpd95_4(1:end,2));
else
    z = 0;
    for i = 1:length(ind2)
        z = z + 1;
        indy2(z,1) = ind2(i);
        z = z + 1;
        indy2(z,1) = ind2(i)+1;
    end 
    indy2 = [ 1 ; indy2; length(hpd95_4(:,1)) ];
    z=0;
    for i = 2:2:length(indy2)
        z = z+1;
        p95_4(z,1) = hpd95_4(indy2(i),1);
        p95_4(z,2) = hpd95_4(indy2(i-1),1);
        p95_4(z,3) = sum(hpd95_4(indy2(i-1):indy2(i),2));    
    end
    p95_4 = flipud(p95_4);
end

% Convert output to BC/BCE if necessary
if strcmpi(yeartype,'BCE/CE') == 1
    prob(:,1) = (prob(:,1)-1950) * -1;
    p95_4(:,1:2) = (p95_4(:,1:2)-1950) * -1;
    p68_2(:,1:2) = (p68_2(:,1:2)-1950) * -1;
end

%%%%% ----- Start of plotting module ----- %%%%%

if strcmpi(yeartype,'Cal BP') == 1
    
    prob2 = prob(cumsum(prob(:,2)) > 0.001 & cumsum(prob(:,2)) < 0.999);
    yrrng = (prob2(end,1) - prob2(1,1))/2;
    
    % round to nearest hundred for nice plot limits
    syr = (10^2) * round((prob2(1,1)-yrrng) / (10^2));
    eyr = (10^2) * round((prob2(end,1)+yrrng) / (10^2));
    
    ind = find(curvecal >= syr & curvecal <= eyr);
    curvecal = curvecal(ind);
    curve14c = curve14c(ind);
    curve1sig = curve1sig(ind);
    
elseif strcmpi(yeartype,'BCE/CE') == 1
    
    prob2 = prob(cumsum(prob(:,2)) > 0.001 & cumsum(prob(:,2)) < 0.999);
    prob2 = flipud(prob2);
    yrrng = (prob2(end,1) - prob2(1,1))/2;
   
    % round to nearest hundred for nice plot limits
    syr = (10^2) * round((prob2(1,1)-yrrng) / (10^2));
    eyr = (10^2) * round((prob2(end,1)+yrrng) / (10^2));
    
    curvecal = (curvecal-1950) * -1;
    ind = find(curvecal >= syr & curvecal <= eyr);
    curvecal = curvecal(ind);
    curve14c = curve14c(ind);
    curve1sig = curve1sig(ind);
      
end

if plotme==1
    figure(14)
    clf

    %----- Plot ProbDistFunc
    axpdf = axes;
    axes(axpdf)
    area(prob(:,1),prob(:,2),'edgecolor','none')
    axpdfylims = ylim;
    axpdfxlims = xlim;
    area(prob(:,1),prob(:,2)*0.2,'edgecolor',[0 0 0],'facecolor',[0.9 0.9 0.9])
    hold on
    [M N] = size(p95_4);
    for i = 1:M
        if strcmpi(yeartype,'Cal BP') == 1
            area( prob(prob(:,1) <= p95_4(i,1) & prob(:,1) >= p95_4(i,2),1)  , prob(prob(:,1) <= p95_4(i,1) & prob(:,1) >= p95_4(i,2),2)*0.2,'edgecolor','none','facecolor',[0.5 0.5 0.6])
        elseif strcmpi(yeartype,'BCE/CE') == 1
            area( prob(prob(:,1) >= p95_4(i,1) & prob(:,1) <= p95_4(i,2),1)  , prob(prob(:,1) >= p95_4(i,1) & prob(:,1) <= p95_4(i,2),2)*0.2,'edgecolor','none','facecolor',[0.5 0.5 0.6])
        end
    end
    %----- Plot cal curve
    axcurve = axes;
    axes(axcurve)
    xdata = curvecal;
    ydata = curve14c;
    onesig = curve1sig;
    fill([xdata' fliplr(xdata')],[ydata'+onesig' fliplr(ydata'-onesig')],[0.6 0.6 0.6],'edgecolor','none');
    hold on
    axcurveylims = ylim;
    axcurvexlims = [min(curvecal) max(curvecal)];

    %----- Plot raw data if intcal13 is selected
    if strcmpi('intcal13',calcurve) == 1;
        
        axraw = axes;
        axes(axraw)
        
        rd = load('private/IntCal13 raw data.txt');
        
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
        
        raw13_cal = rd(ind,3);
        raw13_calsigma = rd(ind,5);
        raw13_14c = rd(ind,6);
        raw13_14csigma = rd(ind,7);
        
        for i = 1:length(raw13_cal)       
           
            % x error bars
            plot([raw13_cal(i)-raw13_calsigma(i) raw13_cal(i)+raw13_calsigma(i)],[raw13_14c(i) raw13_14c(i)],'-','color',[0.8 0.8 0.8])
            
            if i == 1
                hold on
            end
                    
            % y error bars
            plot([raw13_cal(i) raw13_cal(i)],[raw13_14c(i)-raw13_14csigma(i) raw13_14c(i)+raw13_14csigma(i)],'-','color',[0.8 0.8 0.8])
            
            ymaxes(i) = raw13_14c(i)+raw13_14csigma(i);
            ymins(i) = raw13_14c(i)-raw13_14csigma(i);
            
        end
        
        axrawylims = [min(ymins) max(ymaxes)];
        axrawxlims = xlim;
        
    end

    %----- Plot 14C age normal distribution(s)
    axgauss = axes;
    axes(axgauss);

    gaussrange = [c14age-4*c14err:c14age+4*c14err];
    gauss = normpdf(gaussrange,c14age,c14err);

    gaussrangeorig = [c14ageorig-4*c14errorig:c14ageorig+4*c14errorig];
    gaussorig = normpdf(gaussrangeorig, c14ageorig, c14errorig);

    area(gauss, gaussrange);
    axgaussylims = ylim;
    axgaussxlims = xlim;
    area(gauss*0.2, gaussrange,'edgecolor','none','facecolor',[0.7 0.4 0.4]);
    hold on
    plot(gaussorig*0.2, gaussrangeorig,'linestyle','-','color',[0.7 0.4 0.4]);

    %----- set plot settings by axis, starting from back layer to front layer

    axes(axcurve)
    xlim(axcurvexlims)
    if strcmpi(yeartype, 'Cal BP') == 1
        set(gca, 'XDir', 'reverse')
    end
    set(gca,'color','none')
    lab1 = ylabel('^1^4C yr BP');
    lab2 = xlabel(yearlabel);
    set( gca, 'TickDir', 'out' );
    if strcmpi('intcal13',calcurve) == 1;
        ylim(axrawylims)
    else
        ylim(axcurveylims)
    end
    yt=get(gca,'ytick');
    ytl=textscan(sprintf('%1.0f \n',yt),'%s','delimiter','');
    set(gca,'yticklabel',ytl{1})
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
    if strcmpi(yeartype, 'Cal BP') == 1
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
    if strcmpi('intcal13',calcurve) == 1;
        ylim(axrawylims)
    else
        ylim(axcurveylims)
    end

    if strcmpi('intcal13',calcurve) == 1;
        axes(axraw)
        hold on
        if strcmpi(yeartype, 'Cal BP') == 1
            set(gca, 'XDir', 'reverse')
        end
        set(gca,'color','none')
        xlim(axcurvexlims)
        ylim(axrawylims)
        set(gca,'xticklabel',[]);
        set(gca,'xtick',[]);
        set(gca,'yticklabel',[]);
        set(gca,'ytick',[]);
    else
        axes(axpdf)
    end

    %----- Plot some text on the final axis

    if strcmpi(yeartype, 'Cal BP') == 1
        
        xlims = xlim;
        ylims = ylim;
        
        xaxrange = xlims(2) - xlims(1);
        yaxrange = ylims(2) - ylims(1);
        
        [M N] = size(p95_4);
        
        text(xlims(2)-0.63*xaxrange, ylims(2)-0.03*yaxrange, ['^1^4C date: ',num2str(c14ageorig),' \pm ',num2str(c14errorig),' ^1^4C yr BP'])
        
        if extralabel == 1;
            text(xlims(2)-0.63*xaxrange, ylims(2)-0.06*yaxrange, [reslabel,': ',num2str(resage),' \pm ',num2str(reserr), ' ^1^4C yr'])
            ypos = 0.11;
        else
            ypos = 0.08;
        end
        
        if M == 1
            text(xlims(2)-0.63*xaxrange, ylims(2)-ypos*yaxrange, ['Cal age 95.45% HPD interval:'])
        else
            text(xlims(2)-0.63*xaxrange, ylims(2)-ypos*yaxrange, ['Cal age 95.45% HPD intervals:'])
        end
        
        for i = 1:M
            ypos = ypos+0.03;
            text(xlims(2)-0.63*xaxrange, ylims(2)-ypos*yaxrange, [num2str(floor(p95_4(i,3)*1000)/10),'%: ',num2str(p95_4(i,1)),' to ',num2str(p95_4(i,2)),' cal yr BP'])
        end
        
        text(xlims(2)-0.02*xaxrange, ylims(2)-0.03*yaxrange, ['MatCal 1.0 (B.C. Lougheed)'])
        text(xlims(2)-0.02*xaxrange, ylims(2)-0.06*yaxrange, [calcurve, ' ', cite])
        
    elseif strcmpi(yeartype,'BCE/CE')
        
        xlims = xlim;
        ylims = ylim;
        
        xaxrange = abs( xlims(1) - xlims(2) );
        yaxrange = ylims(2) - ylims(1);
        
        [M N] = size(p95_4);
        
        text(xlims(1)+0.63*xaxrange, ylims(2)-0.03*yaxrange, ['^1^4C date: ',num2str(c14ageorig),' \pm ',num2str(c14errorig),' ^1^4C yr BP'])
        
        if extralabel == 1;
            text(xlims(1)+0.63*xaxrange, ylims(2)-0.06*yaxrange, [reslabel,': ',num2str(resage),' \pm ',num2str(reserr), ' ^1^4C yr'])
            ypos = 0.11;
        else
            ypos = 0.08;
        end
        
        if M == 1
            text(xlims(1)+0.63*xaxrange, ylims(2)-ypos*yaxrange, ['Cal age 95.45% HPD interval:'])
        else
            text(xlims(1)+0.63*xaxrange, ylims(2)-ypos*yaxrange, ['Cal age 95.45% HPD intervals:'])
        end
        
        for i = 1:M
            ypos = ypos+0.03;
            text(xlims(1)+0.63*xaxrange, ylims(2)-ypos*yaxrange, [num2str(floor(p95_4(i,3)*1000)/10),'%: ',num2str(p95_4(i,1)),' to ',num2str(p95_4(i,2)),' ',yearlabel])
        end
        
        text(xlims(1)+0.02*xaxrange, ylims(2)-0.03*yaxrange, ['MatCal 1.0 (B.C. Lougheed)'])
        text(xlims(1)+0.02*xaxrange, ylims(2)-0.06*yaxrange, [calcurve, ' ', cite])
        
    end

    if prob(1,2) > 0.000001 || prob(end,2) > 0.000001
        title('Warning! Calibrated age may exceed limits of calibration curve.')
    end

    %----- Fix all font sizes

    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    set(lab1,'FontSize',10)
    set(lab2,'FontSize',10)

    %----- Prep plot for export

    % set paper size (cm)
    set(gcf,'PaperUnits','centimeters')
    Y = 18;
    X = 18;
    set(gcf, 'PaperSize',[X Y])
    % set figure size (cm)
    xSize = 16;
    ySize = 16; 
    % put figure in centre of paper
    xLeft = (X-xSize)/2;
    yBottom = (Y-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yBottom xSize ySize])
    % make background white
    set(gcf,'InvertHardcopy','on');
    set(gcf,'color',1*[1 1 1]);

    % uncomment line below to automatically print Adobe PDF to working directory
    % print(figure(14), '-dpdf', ['MatCal ',num2str(c14ageorig),'±',num2str(c14errorig),'.pdf']);
end

end





