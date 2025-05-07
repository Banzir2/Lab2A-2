function [dV, dR] = Hioki(V, R, ac)
    if ac == 0
        dV = 0.006 * V + 2;
    else
        dV = 0.022 * V + 5;
    end
    if V < 400
        if ac == 0
            dV = 0.006 * V + 0.2;
        else
            dV = 0.02 * V + 0.2;
        end
    end
    if V < 40
        if ac == 0
            dV = 0.006 * V + 0.02;
        else
            dV = 0.02 * V + 0.02;
        end
    end
    if V < 4
        if ac == 0
            dV = 0.006 * V + 0.002;
        else
            dV = 0.02 * V + 0.002;
        end
    end
    if V < 0.4
        if ac == 0
            dV = 0.006 * V + 0.0002;
        else
            dV = 0.02 * V + 0.001;
        end
    end
    
    
    dR = 0.02 * V + 3 * 10000;
    if R < 4e6
        dR = 0.012 * V + 3 * 1000;
    end
    if R < 4e5
        dR = 0.006 * V + 3 * 100;
    end
    if R < 4e4
        dR = 0.006 * V + 3 * 10;
    end
    if R < 4e3
        dR = 0.006 * V + 3;
    end
    if R < 4e2
        dR = 0.006 * V + 3 * 0.1;
    end
end

