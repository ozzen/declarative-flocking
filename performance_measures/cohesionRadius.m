% Olfati-saber
function R = cohesionRadius(q)
R = max(sqrt(sum((q - mean(q,1)).^2,2)));
end