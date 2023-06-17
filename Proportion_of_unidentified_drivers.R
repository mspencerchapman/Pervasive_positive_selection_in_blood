##Infer the proportion of driver mutations that are as yet unidentified

#Function to take a dNdS value
prop_unidentified=function(dNdS,total_nonsynonymous_muts,identified_drivers) {
  implied_total_drivers=total_nonsynonymous_muts*(dNdS-1)/dNdS
  cat(paste("Total implied number of drivers is",round(implied_total_drivers)),sep="\n")
  
  prop_unidentified_=(implied_total_drivers-identified_drivers)/implied_total_drivers
  return(prop_unidentified_)
}

#Apply function to the dNdS values from the data
data_dNdS_with_CI=c(1.11,1.13,1.58)
names(data_dNdS_with_CI)=c("lower_CI","median","upper_CI")
sapply(data_dNdS_with_CI,function(x) prop_unidentified(dNdS=x,
                                                       total_nonsynonymous_muts=39083,
                                                       identified_drivers=1541))
