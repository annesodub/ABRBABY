function [DEV1_se, DEV2_se, STD1_se, STD2_se] = extract_standard_error_DEV_STD(DEV1_avg, DEV2_avg, STD1_avg, STD2_avg)

    %dimensions XXX.data input (for 1 subject) = subjects x channels x timepoints
    %dimensions datasets (.set) output (std) =  channels x timepoints
    
    DEV1_se = std(DEV1_avg,1,1)/sqrt(size(DEV1_avg,1,1));
    
    DEV2_se = std(DEV2_avg,1,1)/sqrt(size(DEV2_avg,1,1));
    
    STD1_se = std(STD1_avg,1,1)/sqrt(size(DEV1_avg,1,1));
    
    STD2_se = std(DEV2_avg,1,1)/sqrt(size(DEV2_avg,1,1));
    
end
