function r_Taylor = rnTaylor(n)
    syms a e theta;
    
    r_num = (a * (1 - e^2))^n;
    r_den = (1 + e * cos(theta))^n;
    taylor_den = taylor(r_den, e, Order = 2);
    r_Taylor = simplify(r_num / taylor_den, 'IgnoreAnalyticConstraints', true);
end