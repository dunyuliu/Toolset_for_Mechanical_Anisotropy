function [mismatch, dipOfStressAxis, dipOfStrainRateAxis] = funcCalcAlignments(stressVec, strainRateVec)
    % function to calculate misalignment and dips of two vectors
    if stressVec(1)<0
        stressVec = - stressVec;
    end
    if strainRateVec(1)<0
        strainRateVec = -strainRateVec;
    end
    mismatch=acosd(dot(stressVec,strainRateVec)/(norm(stressVec)*norm(strainRateVec)));
    dipOfStressAxis = asind(stressVec(3)/norm(stressVec));
    dipOfStrainRateAxis = asind(strainRateVec(3)/norm(strainRateVec));

    if real(mismatch)>90
        mismatch = 180-real(mismatch);
    else
        mismatch = real(mismatch);
    end
    %disp(mismatch)
    %disp(dipOfStressAxis)
    %disp(dipOfStrainRateAxis)
    