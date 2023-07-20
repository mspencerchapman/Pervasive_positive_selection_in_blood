# Pervasive_positive_selection_in_blood
Code accompanying the manuscript 'Pervasive positive selection in blood in 200,618 individuals and novel drivers of clonal haematopoiesis'.

# Running_dnds_on_colony_WGS.R
The code in script 'Running_dnds_on_colony_WGS.R' relates to the section 'Novel FI-driver gene mutations present in haematopoietic stem/progenitor cells" and the generation of the data in Tables S7 and S9.  The plots generated in the folder 'plots/' are not included in the manuscript itself but are graphical representations of the dNdS output the is described in the text.  Note that not all the raw data required to run this code is included here as the main studies generating this data are not all yet published. However it is included such that the methodology can be examined. The output from the dNdS is found in table format in the folder 'data/'.

# Proportion_of_unidentified_drivers.R
The code in script 'Proportion_of_unidentified_drivers.R' relates to the point in the section 'The majority of clonal expansions remain unexplained' regarding the proportion of drivers implied by dNdS that are now recognised when including the novel genes. This is a simple function that shows how the the output values from the global dNdS function, the numbers of identified drivers in the same mutation set.
