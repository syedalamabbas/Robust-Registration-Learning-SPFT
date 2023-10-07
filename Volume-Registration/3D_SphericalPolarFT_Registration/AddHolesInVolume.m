function [ OriginalVolume, OverlapPercentage ] = AddHolesInVolume( maxAddZerosStage, offsetOfVoxels, OriginalVolume )
%ADDHOLESINVOLUME Specialized function randomly add holes in a volume at all layers
% Output holes filled volume
untouchedVolume = OriginalVolume;
[N,~,~] = size(OriginalVolume);
MaxLimitVoxels = (N-1)/2-offsetOfVoxels;   % This is where we add randomly holes by setting certain values in the volume to zero
for iterate =1 :maxAddZerosStage   % For each layer just set the values to zero 
    for s = 1:N
        J = OriginalVolume(:,:,s);
        randomIntegerL = randi(MaxLimitVoxels);
        randomIntegerR = randi(MaxLimitVoxels);
        J((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR,(N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR) = zeros(length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR),length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR));
        OriginalVolume(:,:,s) = J;
    end
end

%% Computing the overlap percentage now
NumberOfValuesInOriginal  = sum(sum(sum(logical(untouchedVolume))));
NumberOfValuesInHolesObject  = sum(sum(sum(logical(OriginalVolume))));
OverlapPercentage =  (NumberOfValuesInHolesObject / NumberOfValuesInOriginal ) * 100;
disp(['Overlap Percentage at stage ', num2str(maxAddZerosStage),'  is ', num2str(OverlapPercentage)]);


end

