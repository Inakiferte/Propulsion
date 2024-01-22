function dY_dx = Func (x,Y,m_e,m_i,v_s,gammad_hat,n_s,gammam_hat,cn_hat,nu_s,sigma_en,sigma_0,w_c,alpha_B,cn,e,Mode)
    %===================================
    % INPUTS:
    %       Y(1) = vi_hat
    %       Y(2) = ne_hat
    %       Y(3) = Te_hat
    %       Y(4) = phi_hat
    %===================================
    
    % Create solution vector
    dY_dx = zeros(4,1);
    
    %========================================================
    % Compute Magnetic field depending the chosen mode
    % Mode = 1 : B = cte.
    % Mode = 2 : B = Gaussian.
    %========================================================

    if Mode==1
        wc = w_c;
    else
        B  = Bx(x);
        wc = (e * B) / m_e;              % [rad/s]
    end
    

    % Compute auxiliary values
    aux   = (5 / 3) * Y(3) - Y(1)^2;
    m_aux = m_e / m_i;
    vex_hat = (Y(1) - (gammad_hat / Y(2)));
    n_n_hat = ((gammam_hat / cn_hat) - (Y(2) * (Y(1) / cn_hat)));
    nuen_hat = ((gammam_hat / cn_hat) - (Y(2) * (Y(1) / cn_hat)) * (sigma_en / sigma_0) * sqrt(8 * Y(3) / pi));
    wc_hat = wc / nu_s;
    nue_hat = wc_hat^2 / (nuen_hat + (alpha_B * wc_hat));
    nuion_hat = n_n_hat * sqrt(8 * Y(3) / pi) * (1 + 2 * Y(3)) * exp(-1 / Y(3));
    Ei_p = 2.5;
    vn_hat = cn / v_s; 


    % Add the function
    dY_dx(1) = aux^-1 * ((5 / 3) * m_aux * (Y(1) - (gammad_hat / Y(2))) * nue_hat * Y(1) ...
        + nuion_hat * ((5 / 3) * Y(3) + Y(1) * (Y(1) - vn_hat) - (Y(1) / vex_hat) * ...
        ((2 / 3) * Ei_p + (5 / 3) * Y(3))));

    dY_dx(2) = aux^-1 * (-(5 / 3) * m_aux * Y(2) * vex_hat * nue_hat - Y(2) * ...
        nuion_hat * ((2 * Y(1) - vn_hat) - ((2 / 3) * Ei_p + (5 / 3) * Y(3)) * (1 / vex_hat)));

    dY_dx(3) = aux^-1 * (-(2 / 3) * Y(1)^2 * m_aux * vex_hat * nue_hat - ...
        nuion_hat * ((2 / 3) * Y(3) * (2 * Y(1) - vn_hat) - ((Y(1)^2 - Y(3)) / vex_hat) * ...
        ((2 / 3) * Ei_p + (5 / 3) * Y(3))));

    dY_dx(4) = aux^-1 * (-(5 / 3) * m_aux * vex_hat * nue_hat * Y(1)^2 - ...
        nuion_hat * ((5 / 3) * Y(3) * (2 * Y(1) - vn_hat) - (Y(1)^2 / vex_hat) * ...
        ((2 / 3) * Ei_p + (5 / 3) * Y(3))));
end