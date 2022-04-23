function [sector4Eff, sector4Res] = sector4Calc(numBiot,numOmega,numAlpha,numGamma,numBeta1,numBeta2)
w = 1-max(numBeta2,numOmega);
mLA=numBeta1*numAlpha*sqrt(numBiot);
sec4AEff = tanh(mLA)/mLA;
hE=(1+(sec4AEff*numBeta1)/(1-numBeta1));

mLB = numAlpha * numGamma*(1-w)*sqrt(numBiot*(hE));
sec4BEff = tanh(mLB)/mLB;
sector4Eff = sec4BEff*(1+(numBeta1*(sec4AEff-1)));
sector4Res = inv(sector4Eff*numBiot*numAlpha^2*numGamma*(1-w));
end