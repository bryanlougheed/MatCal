function pranges = hpdcalc(ages, probs, ageres, siglevel)
% pranges = hpdcalc(ages, probs, ageres, siglevel)
%
% Calculates highest posterior density (HPD) credible interval age.
%
% Input
% =====
% ages     = vector of ages, equally spaced
% probs    = vector of probabilities for each age
%            Total probability should sum to 1.
% ageres   = the age spacing in ages, e.g. 1 for 1 year
% siglevel = Desired significance level, e.g. 0.95 for 95%
%
% Output
% ======
% pranges = n by 3 matrix
%           Age interval(s) calculated using highest posterior density.
%           Each row contains a probability range in Cols 1 and 2, and
%           the associated probability for that range in Col 3.
%           Probabilities are normalised to between zero and one.
%
% ---------------------------------
% B.C. Lougheed 2016, updated 2021.


% check input
if isvector(ages) == 0 || isvector(probs) == 0
	error('Check that ages and probs are both vectors')
end
if numel(ages) ~= numel(probs)
	error('Check that ages and probs have same number of elements')
end

hpd = [reshape(ages,numel(ages),1) reshape(probs,numel(probs),1)];
hpd = sortrows(hpd, 2);
hpd(:,3) = cumsum(hpd(:,2));
hpd = hpd(hpd(:,3) >= 1-siglevel, :);
hpd = sortrows(hpd,1);
ind1 = find(diff(hpd(:,1)) > ageres);
if isempty(ind1)
	pranges(1,1) = hpd(end,1);
	pranges(1,2) = hpd(1,1);
	pranges(1,3) = sum(hpd(1:end,2));
else
	ind2 = NaN(length(ind1)*2,1);
	for i = 1:length(ind1)
		ind2(i*2-1,1) = ind1(i);
		ind2(i*2,1) = ind1(i)+1;
	end
	ind2 = [ 1 ; ind2; length(hpd(:,1)) ];
	pranges = NaN(length(2:2:length(ind2)),3);
	for i = 2:2:length(ind2)
		pranges(i/2,1) = hpd(ind2(i),1);
		pranges(i/2,2) = hpd(ind2(i-1),1);
		pranges(i/2,3) = sum(hpd(ind2(i-1):ind2(i),2));
	end
	pranges = flipud(pranges);
end

end