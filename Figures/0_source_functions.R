label_group <- function(genotype){
  if(genotype %in% c("TG","Case","CASE")){group = "AD"}else{
    if(genotype %in% c("WT","Control","CONTROL")){group = "Control"}}
  return(group)
}

label_colour <- function(genotype){
  if(genotype %in% c("AD","Case","CASE")){colour = wes_palette("Royal1")[1]}else{
    if(genotype %in% c("Control","CONTROL")){colour = wes_palette("Royal1")[2]}}
  return(colour)
}
