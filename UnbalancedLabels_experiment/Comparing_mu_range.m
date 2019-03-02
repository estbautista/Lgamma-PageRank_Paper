function [mcc,idx_mcc] = Comparing_mu_range(Lap,mu,y,gt)
	
Lg = Lap;
mu_length = length(mu);
parfor m = 1 : mu_length
	[mcc_out, ~] = PageRank_comparing(Lg,mu(m),y,gt);
	mcc_tmp(m) = mcc_out;
end

[mcc,idx_mcc] = max(mcc_tmp);