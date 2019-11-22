%define function to calculate Fano factor
function fano = FanoFactor(x,y)
    % x and y are allele-wise expression vectors
    intrinsic_noise_squared = mean((x-y).^2)/(2*(mean(x))*(mean(y)));
    protein_level = mean(x+y);
    fano = intrinsic_noise_squared*protein_level;
end