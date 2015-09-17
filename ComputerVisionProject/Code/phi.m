function val = phi(point, mu, covar)
    val = exp(-1*(point-mu)^2/(2*covar))/sqrt(2*pi*covar);
end