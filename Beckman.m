function [dI] = Beckman(I, ac)
    if ac == 0
        if I < 0.2
            dI = 0.01 * I + (oom(I) / 1000);
        else
            dI = 0.02 * I + 2 * (oom(I) / 1000);
        end
    else
        if I < 0.2
            dI = 0.015 * I + 3 * (oom(I) / 1000);
        else
            dI = 0.025 * I + 4 * (oom(I) / 1000);
        end
    end
end

