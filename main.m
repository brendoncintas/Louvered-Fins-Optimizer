%%
% EML 5152: Intermediate Heat Transfer - Prof. Chengxian "Charlie" Lin
% Brendon Cintas - 6158218
% Final Project
% Florida International University

%% Section 0 clears the script.
clear
%% Section 1 sets up all input variables before use based on the parametric study done in the paper.

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

h1P = "What is the convection coefficient (h1)?\n";
h1 = input(h1P);
bP = "What is the element width (b)?\n";
b = input(bP);
cP = "What is the element height (c)?\n";
c = input(cP);
TP = "What is the constant temperature boundary (°C)?\n";
T = input(TP);

outputP = "Print outputs? Y/N [Y]\n";
output = input(outputP, "s");
if isempty(output)
    output = 'Y';
end
% Derived values
b1 = b*numOmega; % Louver section width
t = c/numAlpha; % Element thickness
c1 = c*numBeta1; % Louver element height
b2 = numBeta2*b;

%%

%Sector 4 Calculation
[sector4Eff, sector4Res] = sector4Calc(numBiot,numOmega,numAlpha,numGamma,numBeta1,numBeta2,h1);
sector3Res = sector3Calc(numBiot,numPhi,numOmega,numAlpha,numGamma,numBeta1,numBeta2,h1,sector4Eff);
sector2Res = sector2Calc(numBiot,numPhi,numOmega,numAlpha,numGamma,numBeta1);
sector1Res = sector1Calc(numBiot,numOmega,numAlpha,numGamma,numBeta1);

% The thermal resistance model can be found from an analysis of the network
% of resistance linking the sectors to one another, and then through the
% sink and the source. This is dependent on whether there are small bonds
% (Omega < B2) or large bonds (Omega > B2)

if numOmega < numBeta2
    totalRes = sector1Res(1) + (1)/((1/sector1Res(2))+(1/(sector1Res(3)+sector3Res))+(1/(sector1Res(1)+sector2Res)));
end
if numOmega > numBeta2
    I = zeros(1,5);
    I(1) = sector1Res(1);
    I(2) = sector3Res(4);
    I(3) = sector3Res(3) + sector1Res(3);
    I(4) = inv((1/sector3Res(1))+(1/(sector3Res(4)+sector3Res(1)))+(1/(sector3Res(3)+sector4Res)));
    I(5) = inv((1/sector1Res(2))+(1/(sector2Res+sector1Res(1))));
    
    totalRes = ((I(1)+I(2))*(I(3)+I(4)+I(5))+ I(3)*(I(4)+I(5)))/(I(1)*I(4)*(I(3)*I(5))+I(1)*I(2)*(I(3)+I(4)+I(5))+I(2)*I(5)*(I(3)+I(4))+I(3)*I(4)*I(5));
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
