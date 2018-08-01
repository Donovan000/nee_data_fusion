%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- LOO Loop -----------------------------------------------------------

for useQflux = [0,1]
    for useBudyko = [0,1]
        for useIGBP = [0,1]
            for catIGBP = [0,1]
                if useIGBP && catIGBP; continue; end
                loo_regressions(useQflux,useIGBP,useBudyko,catIGBP);
            end
        end
    end
end


%% *** END SCRIPT *********************************************************
