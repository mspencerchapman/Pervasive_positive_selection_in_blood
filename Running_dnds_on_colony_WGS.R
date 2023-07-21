library(GenomicRanges)
library(IRanges)
library("Rsamtools")
library("MASS")
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(MutationalPatterns)
library(gridExtra)
library(ape)
options(stringsAsFactors = FALSE)

root_dir="~/R_work/Pervasive_positive_selection_in_blood/"
data_dir<-paste0(root_dir,"data/")
plots_dir<-paste0(root_dir,"plots/")

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

#Source functions needed for the script
resave_plots=F
genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"

#----------------------------------------------------------------
# IMPORT THE DATA
#----------------------------------------------------------------

#Set data file paths
##NOTE THAT THE DATA FOR "TN" (Triple-negative ET), "CML" (Chronic Myeloid Leukemia) and "AML" are not available on git hub currently

myeloid_data_sets=c("FOETAL","CLONAL_HAEM","AGEING","HCT","MPN","TN","CML","AML")

study_files<-sapply(myeloid_data_sets,function(study) grep(study,x=list.files(path = data_dir,pattern=".txt",full.names = T),value=T))

myeloid_studies<-Map(file=study_files,data_set=myeloid_data_sets,function(file,data_set) {
  study_df<-read.delim(file,stringsAsFactors = F,header=T)
  study_df<-study_df%>%
    dplyr::mutate(mut_ref=paste(Chrom,Pos,Ref,Alt,sep="-"),.before="Chrom")%>%
    dplyr::mutate(dataset=data_set)
  return(study_df)
})
names(myeloid_studies)<-myeloid_data_sets

#Add on the PD_no column for datasets that don't have this column
myeloid_studies<-lapply(myeloid_studies,function(df) {
  if(is.null(df$PD_no)){df$PD_no<-stringr::str_split(df$sampleID,pattern = "_",simplify=T)[,1]}
  return(df)
})

#Import the lymphoid data - this is in a slightly different format
lymphoid_SNVs_file<-paste0(data_dir,"LYMPHOID_sites_per_sample.txt")
Ly_muts<-read.delim(lymphoid_SNVs_file)%>%
  dplyr::rename("mut_ref"=ID)%>% #Rename for consistency
  distinct(mut_ref,donor, .keep_all = TRUE)%>% #Only keep muts from the same donor once - these are likely to be early embryonic variants
  mutate(dataset="Lymphoid")

n_colony_df_full<-read.delim(paste0(data_dir,"number_of_colonies_per_individual.txt"))

#----------------------------------------------------------------
# RUNNING dN/dS ACROSS THE DATA
#----------------------------------------------------------------

#These are the gene sets found to be under selection in the UK biobank data analysis

all_ukbb_CH=c("DNMT3A","ASXL1","PPM1D","TET2","ZNF318","SRSF2","TP53","MTA2","YLPM1","ZBTB33",
              "SF3B1","MYD88","SIK3","GNB1","BRCC3","GNAS","SRSF1","CBL","JAK2","DUSP22","CHEK2","ZNF234","SPRED2","IGLL5","BAX",
              "CCDC115","MAGEC3","SH2B3","CCL22","PHIP","KDM6A","SRCAP")
novel_FI_ukbb<-c("SPRED2","MTA2","YLPM1","ZBTB33",
                 "ZNF318","ZNF234","SRSF1",
                 "IGLL5","MYD88","SIK3","CHEK2","MAGEC3",
                 "CCDC115","BAX","SRCAP","SH2B3","CCL22")

##Set up dnds function to run on data in the format we have imported
library(dndscv)
dndscv_on_details=function(details,id="this_sample",outp=1,max_muts_per_gene_per_sample = Inf,use_indel_sites=T,max_coding_muts_per_sample = Inf,...) {
  require(dndscv)
  if(!"sampleID"%in%colnames(details)) {
    details$sampleID=id
  }
  muts=details[,c("sampleID","Chrom","Pos","Ref","Alt")]
  colnames(muts)<-c("sampleID","chr","pos","ref","alt")
  
  dndscvout=dndscv(muts,outp=outp,max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,use_indel_sites=use_indel_sites,min_indels=2,max_coding_muts_per_sample = max_coding_muts_per_sample,...)
  return(dndscvout)
}
mutation_type_vec=c("All","Missense","Nonsense","Splice site","Truncating")
names(mutation_type_vec)=c("wall","wmis","wnon","wspl","wtru")

##Run this across the myeloid data
dndsout_myeloid<-dndscv_on_details(details=dplyr::bind_rows(myeloid_studies)%>%dplyr::mutate(sampleID=PD_no),
                                   outp=3)
dndsout_myeloid_restricted<-dndscv_on_details(details=dplyr::bind_rows(myeloid_studies)%>%dplyr::mutate(sampleID=PD_no),
                                              gene_list=novel_FI_ukbb,
                                              outp=3)
n_colony_df_full%>%
  dplyr::summarise(total_colonies=sum(n_colonies),n_ind=nrow(.))

qval_cutoff=0.1
dndsout_myeloid_restricted$sel_cv%>%
  filter(qmis_cv<qval_cutoff|qtrunc_cv<qval_cutoff|qallsubs_cv<qval_cutoff|qglobal_cv<qval_cutoff)%>%
  arrange(qallsubs_cv)%>%dplyr::select(1:6)

##Run this across the lymphoid data
dndsout_lymphoid<-dndscv_on_details(details=Ly_muts,outp=3)
dndsout_lymphoid_restricted<-dndscv_on_details(details=Ly_muts,gene_list=all_ukbb_CH,outp=3)

dndsout_lymphoid_restricted$sel_cv%>%
  filter(qmis_cv<qval_cutoff|qtrunc_cv<qval_cutoff|qallsubs_cv<qval_cutoff)%>%
  arrange(qglobal_cv)

dndsout_lymphoid_restricted$sel_cv

#----------------------------------------------------------------
# DISPLAY & SAVE SUMMARY TABLES OF THESE RESULTS
#----------------------------------------------------------------

dataset_key<-dplyr::bind_rows(myeloid_studies)%>%
  distinct(PD_no,dataset, .keep_all = F)

dndsout_myeloid_restricted$annotmuts%>%
  left_join(dataset_key,by=c("sampleID"="PD_no"))%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  group_by(dataset,gene,impact)%>%
  summarise(n=n())%>%
  pivot_wider(names_from="impact",values_from = "n")%>%
  replace(is.na(.), 0)%>%
  print(n=30)

dndsout_myeloid_restricted$annotmuts%>%
  left_join(dataset_key,by=c("sampleID"="PD_no"))%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  group_by(gene,impact)%>%
  summarise(n=n())%>%
  pivot_wider(names_from="impact",values_from = "n")%>%
  replace(is.na(.), 0)%>%
  print(n=30)

dndsout_myeloid_restricted$annotmuts%>%
  left_join(dataset_key,by=c("sampleID"="PD_no"))%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  dplyr::select(gene,chr,pos,ref,mut,aachange,ntchange,impact,sampleID,dataset)%>%
  arrange(gene,impact,dataset)%>%
  readr::write_csv(file="myeloid_nonsynonymous_muts_in_novel_genes.csv")

dndsout_lymphoid_restricted$annotmuts%>%
  left_join(dataset_key,by=c("sampleID"="PD_no"))%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  dplyr::select(gene,chr,pos,ref,mut,aachange,ntchange,impact,sampleID,dataset)%>%
  arrange(gene,impact,dataset)%>%
  readr::write_csv(file="lymphoid_nonsynonymous_muts_in_novel_genes.csv")

#-------------------------------------------------------------------------
# GET SUMMARY PLOTS OF THE dN/dS OUTPUT AND Q-VALUES FOR NOVEL UKBB GENES
#-------------------------------------------------------------------------

#Import df of the variant types under selection in by UKBB data
UKBB_variant_types_file<-paste0(data_dir,"UKBB_newgene_variant_type.xlsx")
variant_type_df<-readxl::read_excel(UKBB_variant_types_file)
variant_type_df$Truncating_selection<-ifelse(grepl("Nonsense|Indels",variant_type_df$`Selecting variant type`),T,F)
variant_type_df$Missense_selection<-ifelse(grepl("Missense",variant_type_df$`Selecting variant type`),T,F)

variant_type_df$qtrunc<-sapply(variant_type_df$Genes,function(Gene) {
  dndsout_myeloid_restricted$sel_cv%>%
    filter(gene_name==Gene)%>%
    pull(qtrunc_cv)
})

variant_type_df$qmis<-sapply(variant_type_df$Genes,function(Gene) {
  dndsout_myeloid_restricted$sel_cv%>%
    filter(gene_name==Gene)%>%
    pull(qmis_cv)
})

variant_type_df$wmis<-sapply(variant_type_df$Genes,function(Gene) {
  dndsout_myeloid_restricted$sel_cv%>%
    filter(gene_name==Gene)%>%
    pull(wmis_cv)
})

variant_type_df$wnon<-sapply(variant_type_df$Genes,function(Gene) {
  dndsout_myeloid_restricted$sel_cv%>%
    filter(gene_name==Gene)%>%
    pull(wnon_cv)
})

qvalue_threshold=0.1
p1<-variant_type_df%>%
  filter(Truncating_selection)%>%
  arrange(qtrunc)%>%
  mutate(signif=ifelse(qtrunc<0.01,"**",ifelse(qtrunc<qvalue_threshold,"*","")))%>%
  ggplot(aes(x=forcats::fct_reorder(factor(Genes),qtrunc),y=qtrunc,label=signif,col=qtrunc<qvalue_threshold))+
  geom_point()+
  geom_hline(yintercept=qvalue_threshold,linetype=2)+
  theme_bw()+
  my_theme+
  #geom_text(aes(y=0.11))+
  theme(axis.text.x=element_text(angle=90),legend.position = "none")+
  labs(x="Gene",y="dN/dS q-value\n(truncating mutations)")

p2<-variant_type_df%>%
  filter(Missense_selection)%>%
  arrange(qtrunc)%>%
  mutate(signif=ifelse(qmis<0.01,"**",ifelse(qmis<qval_cutoff,"*","")))%>%
  ggplot(aes(x=forcats::fct_reorder(factor(Genes),qmis),y=qmis,col=qmis<qval_cutoff))+
  geom_point()+
  geom_hline(yintercept=qval_cutoff,linetype=2)+
  theme_bw()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),legend.position = "none")+
  labs(x="Gene",y="dN/dS q-value\n(Missense mutations)")

library(gridExtra)
ggsave(filename = paste0(plots_dir,"dNdS_qvalue_plot.pdf"),gridExtra::arrangeGrob(p1,p2,nrow=1,widths=c(11,9)),width=3,height=2)

p3<-variant_type_df%>%
  mutate(signif=ifelse(qtrunc<0.01,"**",ifelse(qtrunc<qval_cutoff,"*","")))%>%
  ggplot(aes(x=wnon,y=-log10(qtrunc),label=Genes,col=qtrunc<qval_cutoff))+
  geom_point(size=0.7)+
  geom_hline(yintercept=-log10(qval_cutoff),linetype=2)+
  # scale_y_log10(breaks=10^seq(-6,0,1))+
  # scale_y_reverse()+
  scale_x_log10(limits=c(1,100))+
  ggrepel::geom_text_repel(show.legend = FALSE,size=1,max.overlaps = 20)+
  theme_bw()+
  my_theme+
  labs(x="dN/dS value (nonsense variants)",y=expression(-log["10"]*" q-value (all truncating variants)"),col="Significant q-value")+
  theme(legend.position="none")

p4<-variant_type_df%>%
  mutate(signif=ifelse(qmis<0.01,"**",ifelse(qmis<qval_cutoff,"*","")))%>%
  ggplot(aes(x=wmis,y=-log10(qmis),label=Genes,col=qmis<qval_cutoff))+
  geom_point(size=0.7)+
  geom_hline(yintercept=-log10(qval_cutoff),linetype=2)+
  #scale_y_log10(breaks=10^seq(-16,0,2))+
  #scale_y_reverse()+
  scale_x_log10(limits=c(1,20))+
  ggrepel::geom_text_repel(show.legend = FALSE,size=1,max.overlaps = 20)+
  theme_bw()+
  my_theme+
  labs(x="dN/dS value (missense variants)",y=expression(-log["10"]*" q-value (missense variants)"),col="Significant q-value")+
  theme(legend.position="none")

ggsave(filename = paste0(plots_dir,"dNdS_by_qvalue_all_plot.pdf"),gridExtra::arrangeGrob(p3,p4,nrow=1,widths=c(1,1)),width=6,height=2.5)

p5<-variant_type_df%>%
  filter(Truncating_selection)%>%
  mutate(signif=ifelse(qtrunc<0.01,"**",ifelse(qtrunc<qval_cutoff,"*","")))%>%
  ggplot(aes(x=wnon,y=-log10(qtrunc),label=Genes,col=qtrunc<qval_cutoff))+
  geom_point(size=0.7)+
  geom_hline(yintercept=-log10(qval_cutoff),linetype=2)+
  # scale_y_log10(breaks=10^seq(-6,0,1))+
  # scale_y_reverse()+
  scale_x_log10(limits=c(1,100))+
  ggrepel::geom_label_repel(show.legend = FALSE,size=1.5)+
  theme_bw()+
  my_theme+
  labs(x="dN/dS value (nonsense variants)",y=expression(-log["10"]*" q-value (all truncating variants)"),col="Significant q-value")+
  theme(legend.position="none")

p6<-variant_type_df%>%
  filter(Missense_selection)%>%
  mutate(signif=ifelse(qmis<0.01,"**",ifelse(qmis<qval_cutoff,"*","")))%>%
  ggplot(aes(x=wmis,y=-log10(qmis),label=Genes,col=qmis<qval_cutoff))+
  geom_point(size=0.7)+
  geom_hline(yintercept=-log10(qval_cutoff),linetype=2)+
  #scale_y_log10(breaks=10^seq(-16,0,2))+
  #scale_y_reverse()+
  scale_x_log10(limits=c(1,20))+
  ggrepel::geom_label_repel(show.legend = FALSE,size=1.5)+
  theme_bw()+
  my_theme+
  labs(x="dN/dS value (missense variants)",y=expression(-log["10"]*" q-value (missense variants)"),col="Significant q-value")+
  theme(legend.position="none")

ggsave(filename = paste0(plots_dir,"dNdS_by_qvalue_plot.pdf"),gridExtra::arrangeGrob(p5,p6,nrow=1,widths=c(1,1)),width=6,height=2.5)

#-------------------------------------------------------------------------
# LOOK DEEPER AT THE IGLL5 MUTATIONS - what colony phenotypes are they in?
#-------------------------------------------------------------------------

##### Establish what cell types the IGLL mutations are in, and how this compares to the overall cell type ratios
dat<-read.delim(paste0(data_dir,"colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt"))

Ly_muts_with_colony_info<-dndsout_lymphoid_restricted$annotmuts%>%
  mutate(mut_ref=paste(chr,pos,ref,mut,sep="-"))%>%
  left_join(Ly_muts%>%mutate(mut_ref=paste(Chrom,Pos,Ref,Alt,sep="-")),by=c("mut_ref"))%>%
  left_join(dat,by=c("sampleID.y"="colony"))

Ly_muts_with_colony_info%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  dplyr::select(gene,chr,pos,ref,mut,aachange,ntchange,impact,"sampleID"=sampleID.x,"CellType"=Cell.type2,dataset)%>%
  filter(CellType!="HSC")%>%
  arrange(gene,impact,dataset)%>%
  readr::write_csv(file=paste0(data_dir,"lymphoid_nonsynonymous_muts_in_novel_genes.csv"))

Ly_muts_with_colony_info%>%
  filter(gene=="IGLL5")%>%
  group_by(Donor)%>%
  dplyr::summarise(n=n())

Ly_muts_with_colony_info%>%
  filter(gene=="IGLL5" & impact!="Synonymous")%>%
  group_by(Donor)%>%
  dplyr::summarise(n=n())

prop.test(x=c(3,14,0,0,0,1),n=c(21,30,7,10,1,5))

#Compare this to the total numbers of colonies of different cell types
dat%>%group_by(Donor,Cell.type2)%>%dplyr::summarise(n=n())%>%filter(Cell.type2=="Memory B")
dat%>%group_by(Donor,Cell.type2)%>%dplyr::summarise(n=n())%>%filter(Cell.type2=="Naive B")

#-------------------------------------------------------------------------
# HOW MANY GENES HAVE AT LEAST ONE NON-SYNONYMOUS MUTATION IN LYMPHOID COLONIES?
#-------------------------------------------------------------------------

Ly_genes<-Ly_muts_with_colony_info%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%pull(gene)%>%unique()

