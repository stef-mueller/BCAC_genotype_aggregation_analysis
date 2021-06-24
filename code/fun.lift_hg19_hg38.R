# --- load necessary libraries
require(liftOver)
require(rtracklayer)
require(tidyverse)


# necessary function to flatten GRange List elements
#--- function taken from package {biovizBase
flatGRL = function(object,indName = "grl_name"){
  idx <- togroup(PartitioningByWidth(object))
  gr <- stack(object, indName)
  values(gr) <- cbind(values(gr), values(object)[idx, , drop = FALSE])
  gr
}


lift_hg19_hg38 = function(data, chain){
  # objective: liftOver positions from hg19 to hg38
  # input:  data [data.frame]: with columns chr (formatted as "chr[1-22,X,Y]") and column pos_hg19
  #         chain [path]: path to downloaded USCS liftOver chain file 
  # output: input data.frame with additional column pos_hg38
  # requirement: libraries: liftOver, tidyverse, rtracklayer, GRanges
  
  # --- load chain file
  chain_hg19_hg38 = import.chain(chain)
  
  # --- format input to format agreeable with function makeGRangesFromDataFrame
  data_hg19 = data %>% 
    mutate(start=pos_hg19,
           end=pos_hg19)
  
  # --- create GRanges object from input data
  GR_hg19 = makeGRangesFromDataFrame(data_hg19)
  
  # --- liftOver: output will be in GRanges List Format
  GRlist_hg38 = liftOver(x = GR_hg19, chain = chain_hg19_hg38)
  
  # add meta data column to prevent removing positions where liftover not working
  mcols(GRlist_hg38)$pos_hg19 = data$pos_hg19
  
  # --- format GRanges List ot GRange 
  GR_hg38 = flatGRL(GRlist_hg38)
  
  # --- format to data frame
  df_hg38 = data.frame(GR_hg38) %>% 
    dplyr::select("pos_hg38" =start, pos_hg19)
  
  # --- add hg38 positions to input data
  out = data_hg19 %>% 
    left_join(df_hg38, by="pos_hg19") %>% 
    dplyr::select(-start,-end)
  
  return(out)
}

# test_input = data.frame(
#   chr = c("chr5", "chr7", "chrX", "chr1"),
#   pos_hg19 = c(55890630, 4548574, 57545454, 37657872),
#   expected_hg38 =c(56594803, 4508943, 57519021, 37192271)
# )
# 
# test_output = lift_hg19_hg38(data = test_input,
#                              chain = "/home/stef/WORK/UCL/UGI/projects/breast_cancer/publication/data/input/hg19ToHg38.over.chain")

