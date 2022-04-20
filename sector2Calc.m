function sector2Res = sector2Calc(numBiot,numPhi,numOmega,numAlpha,numGamma,numBeta1)
% Sector 2 is the portion of the lover panel directly above the bond. This
% sector is assumed to be a fin of uniform base temperature acting at its
% midpoint and whose tip and edges are insulated.

mL = numBeta1*numAlpha*sqrt(numBiot*numPhi);
sector2Eff = tanh(mL)/mL;
sector2Res = (1)/(sector2Eff*numBiot*numPhi*numAlpha^2*numBeta1*numOmega*numGamma);
end