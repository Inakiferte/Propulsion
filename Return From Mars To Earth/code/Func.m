function dY_dt = Func (t,Y,mu,g,F,Isp)
    % Create solution vector
    dY_dt = zeros(5,1);
    
    % Compute angles
    sa = Y(3) / sqrt(Y(3)^2 + Y(4)^2);
    ca = Y(4) / sqrt(Y(3)^2 + Y(4)^2);

    % Add the function
    dY_dt(1) = Y(3);
    dY_dt(2) = Y(4) / Y(1);
    dY_dt(3) = (-mu / Y(1)^2) + (Y(4)^2 / Y(1)) + (F * sa/ Y(5));
    dY_dt(4) = (-Y(3) * Y(4) / Y(1)) + (F * ca / Y(5));
    dY_dt(5) = -F / (Isp * g);
end