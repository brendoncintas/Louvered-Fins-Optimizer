function sector3Res = sector3Calc(numBiot,numPhi,numOmega,numAlpha,numGamma,numBeta1,numBeta2,h1,sector4Eff)

% In sector 3, there are two cases for either a large bond (where Omega >
% Beta2) and small bonds. An initial check using if states checks to see
% what case and then applies calculations based on that. 

if numOmega > numBeta2
    % If it is a large bond, then zone 3b is assumed to act like an element
    % of uniform temperature acting at the centroid. Heat is sourced by
    % conduction from the borders with sector 1 and the bond.
    
    % Conduction resistance from these boundaries expressed in
    % dimensionless form are shown below:
    
    sector31Res = (numOmega-numBeta2)*numGamma / (2*(1-numBeta1));
    sector3BRes = (1-numBeta1)/(2*(numOmega-numBeta2)*numGamma);
    
    % Some heat is also lost through the sink:
   
    sector3SRes = 1/(numBiot*numAlpha^2*(numOmega-numBeta2)*(1-numBeta1)*numGamma);
    % Then, efficiency is calculated. This will be needed in the overall
    % resistance calculated for zone 3a.
    
    mL = numBeta1*numAlpha*sqrt(numBiot);
    sec3AEff = tanh(mL)/mL;
    
    %Finally, zone 3A resistance is calculated.
    sector3ARes = 1/(sec3AEff*numBiot*numAlpha^2*numGamma*numBeta1*(numOmega-numBeta2));
    
    sector3Res = [sector3ARes, sector3SRes, sector31Res, sector3BRes];
end

if numOmega < numBeta2
    % For small bonds, Zone 3B is not bonded to the source. Heat enters
    % through 3B through the border with sector 1. Thus, it can be modeled
    % the same as 4B as a fin having an effective film coefficient.
    
    mLA=numBeta1*numAlpha*sqrt(numBiot*numPhi);
    sec3AEff = tanh(mLA)/mLA;
    hEquiv = h1*(1+(sec3AEff*numPhi*numBeta1)/(1-numBeta1));

    mL = (numBeta2-numOmega)*numGamma*numAlpha*sqrt(numBiot*(hEquiv/h1));
    ht = h1*((sector4Eff*(1-numBeta2)*numGamma*numAlpha)/(1-numBeta1));
    sector3Res = (1-numBeta1)*numAlpha*sqrt(numBiot*(hEquiv/h1))*(tanh(mL)*(csch(mL)^2)/((h1/ht)*sqrt((hEquiv/h1)*(1/numBiot))+tanh(mL))); 
    sector3Res = inv(sector3Res);
end