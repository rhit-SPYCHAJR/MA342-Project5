function [temps, spins, energies] = simulate1D(pts, iter, t0, k, J, H, rate, randStart, plotRes)
%SIMULATE1D simulates the 1D Ising model
%
%inputs: 
%      pts: number of points to simulate
%     iter: number of iterations to simulate 
%       t0: initial temperature
%        k: temperature constant
%        J: charge interaction parameter
%        H: external field parameter
%     rate: rate of heating
%randStart: whether to randomize the start states (1) or not (0)
%  plotRes: whether to create a movie showing the result
%
% outputs:
%      temps: column vector showing the temperature at each step
%      spins: matrix where each row gives the spins at each step 

temps = zeros(iter+1, 1);
spins = ones(iter+1, pts);
energies = zeros(iter+1, 1);
T = t0;
temps(1) = T;

%  randomize the initial spins
if (randStart == 1)
    for c = 1:pts
        spins(1,c) = rSpin();
    end
end

energies(1) = stateEnergy(spins(1,:), J, H);

% simulate states
for i = 2:iter+1
    for j = 1:pts
        reversed = reverse(spins(i-1,:),j);
        currentE = stateEnergy(spins(i-1,:), J, H);
        revE = stateEnergy(reversed, J, H);
        if (revE < currentE)
            spins(i,j) = flip(spins(i-1, j));
        elseif (revE == currentE)
            spins(i,j) = spins(i-1,j);
        else
            diff = currentE - revE;
            prob = exp(diff/T);
            if (rand < prob) 
                spins(i,j) = flip(spins(i-1, j));
            else
                spins(i,j) = spins(i-1,j);
            end
        end
    end
    T = T + rate;
    temps(i) = T;
    energies(i) = stateEnergy(spins(i,:), J, H);
end

%plot total energy
if (plotRes == 1)
    xVals = linspace(0,iter, iter+1);
    cells = num2cell(spins,2);

    energyVals = cellfun(@energy, cells);
    absMagVals = cellfun(@absMag, cells);
    figure(1)
    plot(xVals, energyVals);
    xlabel('Iteration')
    ylabel('Sample energy per point')
    %plot(xVals, heatCapVals);
    figure(2)
    plot(xVals, absMagVals);
    xlabel('Iteration')
    ylabel('Sample abs. mag. per point')
    hold off;
    %plot(xVals, susVals);
end

function [x] = rSpin()
    if (rand < 0.5)
        x = -1;
    else
        x = 1;
    end
end


function [y] = flip(x)
    if x == -1
        y = 1;
    else 
        y = -1;
    end
end

function [energy] = stateEnergy(state, j, h)
    tot = state(1)*state(pts);
    for k = 1:pts-1
        tot = tot + state(k)*state(k+1);
    end
    chargeSum = sum(state);
    energy = -1*j*tot - h*chargeSum;
end

function [x] = energy(state)
    E_T = stateEnergy(state, J, H)/pts;
    x = E_T/pts;
end

function [x] = absMag(state)
    M_bar = sum(state)/pts;
    x = abs(M_bar)/pts;
end

function [result] = reverse(state, pos)
    result = state;
    result(pos) = flip(result(pos));
end



end
