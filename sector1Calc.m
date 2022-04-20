function sector1Res = sector1Calc(numBiot,numOmega,numAlpha,numGamma,numBeta1)

% Sector 1 is the unlouvered fin above the bond. This sector is assumed to
% be an element of uniform temperature acting at its sentroid. Heat is
% conducted into the element from its contacting boundary with the source
% which is at a uniform temperature through resistance

sector1BRes = (1-numBeta1)/(2*numOmega*numGamma);

% Also some of the heat is convected directly to the sink through
% resistance

sector1SRes = 1/(numBiot*numAlpha^2*numOmega*(1-numBeta1)*numGamma);

% Some heat is conducted to the border with sector 2 through resistance R1b
% and to the border with sector 3 through resistance.

sector13Res = (numGamma*numOmega)/(2*numBiot);

% The total sector resistance:

sector1Res = [sector1BRes, sector1SRes, sector13Res];
end