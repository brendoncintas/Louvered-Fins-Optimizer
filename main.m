%%
% EML 5152: Intermediate Heat Transfer - Prof. Chengxian "Charlie" Lin
% Brendon Cintas - 6158218
% Final Project
% Florida International University

%% Section 0 clears the script.
clear
%% Section 1 sets up all input variables before use based on the parametric study done in the paper.
customValuesP = "Use custom values? Y/N [N]";
customValues = input(customValuesP, "s");
if isempty(customValues)
    customValues = 'N';
end

outputP = "Print outputs? Y/N [Y]\n";
output = input(outputP, "s");
if isempty(output)
    output = 'Y';
end

if customValues == 'Y'
    % Parametric values
    numBiotP = "What is the parametric Biot number (h1*t/k)?\n";
    numBiot = input(numBiotP);
    numPhiP = "What is the parametric Phi value (h2/h1)?\n";
    numPhi = input(numPhiP);
    numOmegaP = "What is the parametric Omega value (b1/b)?\n";
    numOmega = input(numOmegaP);
    numAlphaP = "What is the parametric Alpha value (c/t)?\n";
    numAlpha = input(numAlphaP);
    numGammaP = "What is the parametric Gamma value (b/c)?\n";
    numGamma = input(numGammaP);
    numBeta1P = "What is the parametric Beta1 value (c1/c)?\n";
    numBeta1 = input(numBeta1P);
    numBeta2P = "What is the parametric Beta2 value (b2/b)?\n";
    numBeta2 = input (numBeta2P);

    bP = "What is the element width (b)?\n";
    b = input(bP);
    cP = "What is the element height (c)?\n";
    c = input(cP);

end
if customValues == 'N'
    numBiot = 1e-5;
    % numBiot = [1e-4, 1e-5, 1e-6];
    numPhi = 0.6;
    % numPhi = [3.3, 2, 1.6];
    numOmega = 0.5;
    % numOmega = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9];
    numAlpha = 200;
    % numAlpha = [200, 140, 80];
    numGamma = 1.5;
    % numGamma = [1.5, 1.0, 0.5];
    numBeta1 = 0.75;
    %numBeta1 = [0.5:0.01:0.95];
    %numBeta2 = 0.25;
    numBeta2 = [0.25:0.01:0.75];
    b = 10;
    c = 5;
    h1 = 100;
    T = 300;

end
% Derived values
b1 = b*numOmega; % Louver section width
t = c/numAlpha; % Element thickness
c1 = c*numBeta1; % Louver element height
b2 = numBeta2*b;


%%

for i = 1:length(numBeta2)
    % The thermal resistance model can be found from an analysis of the network
    % of resistance linking the sectors to one another, and then through the
    % sink and the source. This is dependent on whether there are small bonds
    % (Omega < B2) or large bonds (Omega > B2)

    if numOmega < numBeta2(i)

        [sector4Eff(i), sector4Res(i)] = sector4Calc(numBiot,numOmega,numAlpha,numGamma,numBeta1,numBeta2(i));
        [sector3ARes(i), sector3SRes(i), sector31Res(i), sector3BRes(i), sector3Res(i)] = sector3Calc(numBiot,numPhi,numOmega,numAlpha,numGamma,numBeta1,numBeta2(i),sector4Eff(i));
        sector2Res(i) = sector2Calc(numBiot,numPhi,numOmega,numAlpha,numGamma,numBeta1);
        [sector1BRes(i), sector1SRes(i), sector13Res(i)] = sector1Calc(numBiot,numOmega,numAlpha,numGamma,numBeta1);
        
        J = inv(sector1SRes(i))+inv(sector13Res(i)+sector3Res(i))+inv(sector1BRes(i)+sector2Res(i));
        totalRes(i) = sector1BRes(i) + inv(J);
    end

    if numOmega > numBeta2(i)

        [sector4Eff(i), sector4Res(i)] = sector4Calc(numBiot,numOmega,numAlpha,numGamma,numBeta1,numBeta2(i));
        [sector3ARes(i), sector3SRes(i), sector31Res(i), sector3BRes(i), sector3Res(i)] = sector3Calc(numBiot,numPhi,numOmega,numAlpha,numGamma,numBeta1,numBeta2(i),sector4Eff(i));
        sector2Res(i) = sector2Calc(numBiot,numPhi,numOmega,numAlpha,numGamma,numBeta1);
        [sector1BRes(i), sector1SRes(i), sector13Res(i)] = sector1Calc(numBiot,numOmega,numAlpha,numGamma,numBeta1);

        I = zeros(1,5);
        I(1) = sector1BRes(i);
        I(2) = sector3BRes(i);
        I(3) = sector31Res(i) + sector13Res(i);
        I(4) = inv(inv(sector3SRes(i))+(inv(sector3BRes(i)+sector3ARes(i)))+(inv(sector31Res(i)+sector4Res(i))));
        I(5) = inv(inv(sector1SRes(i))+(inv(sector2Res(i)+sector1BRes(i))));

        A = (I(1)+I(2))*(I(3)+I(4)+I(5));
        B = I(3)*(I(4)+I(5));
        C = (I(1)*I(4))*(I(3)+I(5));
        D = (I(1)*I(2))*(I(3)+I(4)+I(5));
        E = (I(2)*I(5))*(I(3)+I(4));
        F = I(3)*I(4)*I(5);

        totalRes(i) = (A + B) / (C + D + E + F);
    end
end
if output == 'Y'
    fprintf("Based on input values, your:\n");
    fprintf("Louver section width (b2) is %2f .\n", b2);
    fprintf("Louver section height (c1) is %2f .\n", c1);
    fprintf("Louver section thickness (t) is %2f .\n", 2*t);
    fprintf("Temperature boundary width (b1) is %2f .\n\n", b1);
    fprintf("Total resistance with louvered section is %2f .\n,", totalRes);
end
if output == 'N'
end

hold on
plot(numBeta2,totalRes);

