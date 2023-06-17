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

#my_working_directory<-"/lustre/scratch119/casm/team154pc/ms56/lesion_segregation"
root_dir="~/R_work/Pervasive_positive_selection_in_blood/"


my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch119/realdata/mdt1/team154/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
resave_plots=F

lesion_seg_input_dir="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data"
genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"

#Set data file paths
data_sets=c("MSC_fetal","MF","EM","MSC_BMT")
details_all_list=lapply(data_sets,function(data_set) {
  cat(data_set,sep="\n")
  samples<-readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  details_list<-lapply(samples,function(SampleID) {
    cat(SampleID,sep="\n")
    info<-get_file_paths_and_project(dataset=data_set,Sample_ID=SampleID)
    load(info$filtered_muts_path)
    return(filtered_muts$COMB_mats.tree.build$mat)
  })
  names(details_list)<-samples
  return(details_list)
})
names(details_all_list)<-data_sets

trees_all_list=lapply(data_sets,function(data_set) {
  cat(data_set,sep="\n")
  samples<-readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  trees_list<-lapply(samples,function(SampleID) {
    cat(SampleID,sep="\n")
    info<-get_file_paths_and_project(dataset=data_set,Sample_ID=SampleID)
    tree<-read.tree(info$tree_file_path)
    return(tree)
  })
  names(trees_list)<-samples
  return(trees_list)
})
names(trees_all_list)<-data_sets

#Get number of colonies per donor from the trees
ncolony_df<-Map(dataset_trees=trees_all_list,dataset=names(trees_all_list),function(dataset_trees,dataset) {
  Map(tree=dataset_trees,donor=names(dataset_trees),function(tree,donor) {
    data.frame(dataset=dataset,donor=donor,n_colonies=length(tree$tip.label[!tree$tip.label=="Ancestral"]))
  })%>%dplyr::bind_rows()
})%>%dplyr::bind_rows()

# #Import the MPN data - this is in a slightly different format
# MPN_RDS_files=list.files("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/NW",pattern = ".RDS",full.names = T)
# MPN_samples=gsub(".RDS","",list.files("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/NW",pattern = ".RDS",full.names = F))
# MPN_data<-lapply(MPN_RDS_files,function(RDS_file) {
#   cat(RDS_file)
#   data=readRDS(RDS_file)
#   return(data)
# })
# names(MPN_data)<-MPN_samples
# MPN_details<-Map(sampleID=MPN_samples,data=MPN_data,function(sampleID,data) {
#   temp=data$details%>%
#     #filter(FILTER==0)%>%
#     dplyr::rename("Gene"=GENE)%>%
#     dplyr::select(Chrom,Pos,Ref,Alt,Gene,HGVS_PROTEIN,VC)%>%
#     dplyr::mutate(mut_ref=paste(Chrom,Pos,Ref,Alt,sep="-"),.before=1)%>%
#     dplyr::mutate(sampleID=sampleID)
#   return(temp)
# })
# ncolony_MPN=unlist(lapply(MPN_data,function(data) return(length(data$meta$tree$tip.label)-1)))
# ncolony_MPN_df=data.frame(dataset="MPN",donor=MPN_samples,n_colonies=ncolony_MPN)

#Now import datasets that don't have trees
#Import the lymphoid data - a different format
lymphoid_SNVs_file<-"/lustre/scratch119/casm/team154pc/ms56/dNdS_UKBB/som2_indel_SNV_allDF_v2.txt"
Ly_muts<-read.delim(lymphoid_SNVs_file)%>%
  dplyr::rename("mut_ref"=ID,"Chrom"=chr,"Pos"=pos,"Ref"=ref,"Alt"=mut,"sampleID"=donor,"individual_sample"=sample)%>% #Rename for consistency
  distinct(mut_ref,sampleID, .keep_all = TRUE)%>% #Only keep muts from the same donor once - these are likely to be early embryonic variants
  mutate(dataset="Lymphoid")

#Import the triple-negative ET, SDS and CML data - these all have the same format as each other
studies<-c("MPN","TN","CML","SDS")
root_dir<-"/lustre/scratch119/casm/team154pc/ms56/dNdS_UKBB"
study_files<-sapply(studies,function(study) grep(study,x=list.files(path = root_dir,pattern=".txt",full.names = T),value=T))

other_studies<-Map(file=study_files,data_set=studies,function(file,data_set) {
  study_df<-read.delim(file,stringsAsFactors = F,header=T)
  study_df<-dplyr::select(study_df,-sampleID)%>%
    dplyr::rename("Chrom"=chr,"Pos"=pos,"Ref"=ref,"Alt"=mut,"sampleID"=donor,"Gene"=GENE)%>% #Rename for consistency
    dplyr::mutate(mut_ref=paste(Chrom,Pos,Ref,Alt,sep="-"),.before="Chrom")%>%
    dplyr::mutate(dataset=data_set)
  return(study_df)
})
names(other_studies)<-studies

#catalogue the number of colonies
ncolony_df2<-Map(df=other_studies,dataset=names(other_studies),function(df,dataset){
  df%>%filter(nchild==1)%>%
    group_by(sampleID,node)%>%
    summarise(n=n())%>%
    group_by(sampleID)%>%
    summarise(n_colonies=n())%>%
    mutate(dataset=dataset,.before=1)%>%
    dplyr::rename("donor"=sampleID)
})%>%dplyr::bind_rows()

#n_colony_df_full<-bind_rows(ncolony_df,ncolony_MPN_df,ncolony_df2,data.frame(dataset="AML",donor="tAML",n_colonies=99))
n_colony_df_full<-bind_rows(ncolony_df,ncolony_df2,data.frame(dataset="AML",donor="tAML",n_colonies=99))

#Import the tAML data
tAML_file<-"/lustre/scratch119/casm/team154pc/ms56/dNdS_UKBB/filtered_muts_tAML2_m40_postMS_reduced_a_j_vaf"
load(tAML_file)
tAML_details<-filtered_muts$COMB_mats.tree.build$mat
tAML_details$sampleID<-"tAML"
tAML_details$dataset<-"tAML"

##UKBB analysis
all_ukbb_CH=c("DNMT3A","ASXL1","PPM1D","TET2","ZNF318","SRSF2","TP53","MTA2","YLPM1","ZBTB33",
              "SF3B1","MYD88","SIK3","GNB1","BRCC3","GNAS","SRSF1","CBL","JAK2","DUSP22","CHEK2","ZNF234","SPRED2","IGLL5","BAX",
              "CCDC115","MAGEC3","SH2B3","CCL22","PHIP","KDM6A","SRCAP")
novel_FI_ukbb<-c("SPRED2","MTA2","YLPM1","ZBTB33",
                 "ZNF318","ZNF234","SRSF1",
                 "IGLL5","MYD88","SIK3","CHEK2","MAGEC3",
                 "CCDC115","BAX","SRCAP","SH2B3","CCL22")
##Set up dnds
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

##Run dnds on combined dataset
#combine each dataset
combined_by_dataset<-Map(data_set_details=details_all_list,data_set=names(details_all_list),f=function(data_set_details,data_set) {
  cat(data_set,sep="\n")
  samples<-names(data_set_details)
  dataset_details_combined<-dplyr::bind_rows(Map(details=data_set_details,SampleID=samples,function(details,SampleID) {cat(SampleID,sep="\n");return(details%>%mutate(sampleID=SampleID,dataset=data_set))}))
  return(dataset_details_combined)
})
names(combined_by_dataset)<-names(details_all_list)
combined_by_dataset$Chemo<-combined_by_dataset$EM%>%filter(sampleID%in%c("PX001_2_01","PX002_2_01"))%>%mutate(dataset="Chemo")
combined_by_dataset$EM<-combined_by_dataset$EM%>%filter(!sampleID%in%c("PX001_2_01","PX002_2_01"))%>%mutate(dataset="Normal")

#all_combined_by_dataset<-c(combined_by_dataset,other_studies,list(Ly=Ly_muts,AML=tAML_details,MPN=dplyr::bind_rows(MPN_details)%>%mutate(dataset="MPN")))
all_combined_by_dataset<-c(combined_by_dataset,other_studies,list(Ly=Ly_muts,AML=tAML_details))

##Now combine into single dataset
combined_all<-dplyr::bind_rows(all_combined_by_dataset)
#dndsout_combined<-dndscv_on_details(details=combined_all,outp=3)

###
adult_normal_germline<-c("EM","MPN","MSC_BMT","MSC_fetal","TN","CML","MF","AML","Chemo")
adult_normal_germline_to_include<-c("EM","MPN","MSC_BMT","MSC_fetal","TN","CML","MF","AML")
dndsout_restricted<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[adult_normal_germline_to_include]),
                                      gene_list=novel_FI_ukbb,
                                      outp=3)
n_colony_df_full%>%filter(dataset%in%adult_normal_germline_to_include &
                            !donor%in%c("PX001_2_01","PX002_2_01"))%>%
  dplyr::summarise(total_colonies=sum(n_colonies),n_ind=nrow(.))
qval_cutoff=0.1
dndsout_restricted$sel_cv%>%
  filter(qmis_cv<qval_cutoff|qtrunc_cv<qval_cutoff|qallsubs_cv<qval_cutoff|qglobal_cv<qval_cutoff)%>%
  arrange(qallsubs_cv)%>%dplyr::select(1:6)


###
chemo_exposed<-c("Chemo","AML")
dndsout_chemo_restricted<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[chemo_exposed]),
                                            gene_list=novel_FI_ukbb,
                                            outp=3)

dndsout_chemo_restricted$sel_cv%>%
  filter(qmis_cv<qval_cutoff|qtrunc_cv<qval_cutoff|qallsubs_cv<qval_cutoff)%>%
  arrange(qglobal_cv)

###
lymphoid<-c("Ly")
dndsout_lymphoid<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[lymphoid]),
                                               outp=3)
dndsout_lymphoid_restricted<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[lymphoid]),
                                            gene_list=all_ukbb_CH,
                                            outp=3)

dndsout_lymphoid_restricted$sel_cv%>%
  filter(qmis_cv<qval_cutoff|qtrunc_cv<qval_cutoff|qallsubs_cv<qval_cutoff)%>%
  arrange(qglobal_cv)

dndsout_lymphoid_restricted$sel_cv
##
dataset_key<-combined_all%>%
  distinct(sampleID,dataset, .keep_all = F)

dndsout_restricted$annotmuts%>%
  left_join(dataset_key,by="sampleID")%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  group_by(dataset,gene,impact)%>%
  summarise(n=n())%>%
  pivot_wider(names_from="impact",values_from = "n")%>%
  replace(is.na(.), 0)%>%
  print(n=30)

dndsout_restricted$annotmuts%>%
  left_join(dataset_key,by="sampleID")%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  group_by(gene,impact)%>%
  summarise(n=n())%>%
  pivot_wider(names_from="impact",values_from = "n")%>%
  replace(is.na(.), 0)%>%
  print(n=30)

dndsout_restricted$annotmuts%>%
  left_join(dataset_key,by="sampleID")%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  dplyr::select(gene,chr,pos,ref,mut,aachange,ntchange,impact,sampleID,dataset)%>%
  arrange(gene,impact,dataset)%>%
  readr::write_csv(file="myeloid_nonsynonymous_muts_in_novel_genes.csv")

# dndsout_lymphoid_restricted$annotmuts%>%
#   left_join(dataset_key,by="sampleID")%>%
#   filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
#   dplyr::select(gene,chr,pos,ref,mut,aachange,ntchange,impact,sampleID,dataset)%>%
#   arrange(gene,impact,dataset)%>%
#   readr::write_csv(file="lymphoid_nonsynonymous_muts_in_novel_genes.csv")

#Import df of the variant types under selection in by UKBB data
UKBB_variant_types_file<-"/lustre/scratch119/casm/team154pc/ms56/dNdS_UKBB/UKBB_newgene_variant_type.xlsx"
variant_type_df<-readxl::read_excel(UKBB_variant_types_file)
variant_type_df$Truncating_selection<-ifelse(grepl("Nonsense|Indels",variant_type_df$`Selecting variant type`),T,F)
variant_type_df$Missense_selection<-ifelse(grepl("Missense",variant_type_df$`Selecting variant type`),T,F)

variant_type_df$qtrunc<-sapply(variant_type_df$Genes,function(Gene) {
  dndsout_restricted$sel_cv%>%
    filter(gene_name==Gene)%>%
    pull(qtrunc_cv)
})

variant_type_df$qmis<-sapply(variant_type_df$Genes,function(Gene) {
  dndsout_restricted$sel_cv%>%
    filter(gene_name==Gene)%>%
    pull(qmis_cv)
})

variant_type_df$wmis<-sapply(variant_type_df$Genes,function(Gene) {
  dndsout_restricted$sel_cv%>%
    filter(gene_name==Gene)%>%
    pull(wmis_cv)
})

variant_type_df$wnon<-sapply(variant_type_df$Genes,function(Gene) {
  dndsout_restricted$sel_cv%>%
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
  mutate(signif=ifelse(qmis<0.01,"**",ifelse(qmis<qvalue_threshold,"*","")))%>%
  ggplot(aes(x=forcats::fct_reorder(factor(Genes),qmis),y=qmis,col=qmis<qvalue_threshold))+
  geom_point()+
  geom_hline(yintercept=qvalue_threshold,linetype=2)+
  theme_bw()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),legend.position = "none")+
  labs(x="Gene",y="dN/dS q-value\n(Missense mutations)")

library(gridExtra)
ggsave(filename = "/lustre/scratch119/casm/team154pc/ms56/dNdS_UKBB/dNdS_qvalue_plot.pdf",arrangeGrob(p1,p2,nrow=1,widths=c(11,9)),width=3,height=2)

p3<-variant_type_df%>%
  mutate(signif=ifelse(qtrunc<0.01,"**",ifelse(qtrunc<qvalue_threshold,"*","")))%>%
  ggplot(aes(x=wnon,y=-log10(qtrunc),label=Genes,col=qtrunc<qvalue_threshold))+
  geom_point(size=0.7)+
  geom_hline(yintercept=-log10(qvalue_threshold),linetype=2)+
  # scale_y_log10(breaks=10^seq(-6,0,1))+
  # scale_y_reverse()+
  scale_x_log10(limits=c(1,100))+
  ggrepel::geom_text_repel(show.legend = FALSE,size=1,max.overlaps = 20)+
  theme_bw()+
  my_theme+
  labs(x="dN/dS value (nonsense variants)",y=expression(-log["10"]*" q-value (all truncating variants)"),col="Significant q-value")+
  theme(legend.position="none")

p4<-variant_type_df%>%
  mutate(signif=ifelse(qmis<0.01,"**",ifelse(qmis<qvalue_threshold,"*","")))%>%
  ggplot(aes(x=wmis,y=-log10(qmis),label=Genes,col=qmis<qvalue_threshold))+
  geom_point(size=0.7)+
  geom_hline(yintercept=-log10(qvalue_threshold),linetype=2)+
  #scale_y_log10(breaks=10^seq(-16,0,2))+
  #scale_y_reverse()+
  scale_x_log10(limits=c(1,20))+
  ggrepel::geom_text_repel(show.legend = FALSE,size=1,max.overlaps = 20)+
  theme_bw()+
  my_theme+
  labs(x="dN/dS value (missense variants)",y=expression(-log["10"]*" q-value (missense variants)"),col="Significant q-value")+
  theme(legend.position="none")

ggsave(filename = "/lustre/scratch119/casm/team154pc/ms56/dNdS_UKBB/dNdS_by_qvalue_all_plot.pdf",arrangeGrob(p3,p4,nrow=1,widths=c(1,1)),width=6,height=2.5)

p5<-variant_type_df%>%
  filter(Truncating_selection)%>%
  mutate(signif=ifelse(qtrunc<0.01,"**",ifelse(qtrunc<qvalue_threshold,"*","")))%>%
  ggplot(aes(x=wnon,y=-log10(qtrunc),label=Genes,col=qtrunc<qvalue_threshold))+
  geom_point(size=0.7)+
  geom_hline(yintercept=-log10(qvalue_threshold),linetype=2)+
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
  mutate(signif=ifelse(qmis<0.01,"**",ifelse(qmis<qvalue_threshold,"*","")))%>%
  ggplot(aes(x=wmis,y=-log10(qmis),label=Genes,col=qmis<qvalue_threshold))+
  geom_point(size=0.7)+
  geom_hline(yintercept=-log10(qvalue_threshold),linetype=2)+
  #scale_y_log10(breaks=10^seq(-16,0,2))+
  #scale_y_reverse()+
  scale_x_log10(limits=c(1,20))+
  ggrepel::geom_label_repel(show.legend = FALSE,size=1.5)+
  theme_bw()+
  my_theme+
  labs(x="dN/dS value (missense variants)",y=expression(-log["10"]*" q-value (missense variants)"),col="Significant q-value")+
  theme(legend.position="none")

ggsave(filename = "/lustre/scratch119/casm/team154pc/ms56/dNdS_UKBB/dNdS_by_qvalue_plot.pdf",arrangeGrob(p5,p6,nrow=1,widths=c(1,1)),width=6,height=2.5)


##### Establish what cell types the IGLL mutations are in, and how this compares to the overall cell type ratios
dat<-read.delim("../dNdS_UKBB/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt")

Ly_muts_with_colony_info<-dndsout_lymphoid_restricted$annotmuts%>%
  mutate(mut_ref=paste(chr,pos,ref,mut,sep="-"))%>%
  left_join(dplyr::bind_rows(all_combined_by_dataset[lymphoid])%>%
              mutate(mut_ref=paste(Chrom,Pos,Ref,Alt,sep="-")),by=c("mut_ref"))%>%
  left_join(dat,by=c("individual_sample"="colony"))

Ly_muts_with_colony_info%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  dplyr::select(gene,chr,pos,ref,mut,aachange,ntchange,impact,"sampleID"=sampleID.x,"CellType"=Cell.type2,dataset)%>%
  filter(CellType!="HSC")%>%
  arrange(gene,impact,dataset)%>%
  readr::write_csv(file="lymphoid_nonsynonymous_muts_in_novel_genes.csv")

Ly_muts_with_colony_info%>%
  filter(gene=="IGLL5")%>%
  group_by(Donor)%>%
  dplyr::summarise(n=n())

Ly_muts_with_colony_info%>%
  filter(gene=="IGLL5" & impact!="Synonymous")%>%
  group_by(individual_sample)%>%
  dplyr::summarise(n=n())

prop.test(x=c(3,14,0,0,0,1),n=c(21,30,7,10,1,5))

dat%>%group_by(Donor,Cell.type2)%>%dplyr::summarise(n=n())%>%filter(Cell.type2=="Memory B")
dat%>%group_by(Donor,Cell.type2)%>%dplyr::summarise(n=n())%>%filter(Cell.type2=="Naive B")

########

Ly_genes<-Ly_muts_with_colony_info%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%pull(gene)%>%unique()

dndsout_restricted$annot

########
genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
plot_96profile=function(df,colnames=c("sampleID","chr","pos","ref","alt"),genomeFile) {
  require(Rsamtools)
  require(GenomicRanges)
  require(IRanges)
  mutations = data.frame(sampleID=df[[colnames[1]]],
                         chr=df[[colnames[2]]],
                         pos=df[[colnames[3]]],
                         ref=df[[colnames[4]]],
                         trinuc_ref= "-",
                         mut=df[[colnames[5]]])
  mutations_GRange<-GRanges(mutations$chr, IRanges::IRanges(mutations$pos-1, mutations$pos+1))
  mutations$trinuc_ref = as.vector(Rsamtools::scanFa(file=genomeFile,mutations_GRange ))
  
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  #pdf("Mut_Sig_fetal_shared.pdf")
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full; maxy = max(y)
  h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  }
}
plot_96profile(df=dndsout_lymphoid$annotmuts%>%filter(gene=="IGLL5"),
               colnames=c("sampleID","chr","pos","ref","mut"),
               genomeFile=genome_file
)
###################----------------OLD ANALYSES-------------###################
all_samples<-unique(combined_all$sampleID)
not_published<-c("PX001_2_01","PX002_2_01")
published_samples<-all_samples[!all_samples%in%not_published]
dndsout_combined$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)

#Now do on restricted gene list
dndsout_combined_all_restricted<-dndscv_on_details(details=combined_all,#%>%filter(sampleID%in%published_samples),
                                               gene_list=novel_FI_ukbb,
                                               outp=3)
dndsout_combined_all_restricted$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)

dndsout_combined_all_restricted$sel_cv%>%
  dplyr::select(gene_name,qmis_cv,qtrunc_cv,qallsubs_cv,qglobal_cv)%>%
  dplyr::mutate(signif=qmis_cv<0.1|qtrunc_cv<0.1)%>%
  dplyr::mutate(qmis_cv=sapply(qmis_cv,function(x) max(1e-4,x)),qtrunc_cv=sapply(qtrunc_cv,function(x) max(1e-4,x)))%>%
  ggplot(aes(x=qmis_cv,y=qtrunc_cv,col=signif,label=gene_name))+
  geom_point()+
  scale_x_log10(breaks=10^(-4:0))+
  scale_y_log10(breaks=10^(-4:0))+
  ggrepel::geom_label_repel()+
  theme_classic()+
  geom_hline(yintercept=0.1,linetype=2)+
  geom_vline(xintercept=0.1,linetype=2)
  
#Repeat on published data only
dndsout_combined_published_restricted<-dndscv_on_details(details=combined_all%>%filter(sampleID%in%published_samples),
                                                   gene_list=novel_FI_ukbb,
                                                   outp=3)
dndsout_combined_published_restricted$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)

dndsout_combined_published_restricted$sel_cv%>%
  dplyr::select(gene_name,qmis_cv,qtrunc_cv,qallsubs_cv,qglobal_cv)%>%
  dplyr::mutate(signif=qmis_cv<0.1|qtrunc_cv<0.1)%>%
  dplyr::mutate(qmis_cv=sapply(qmis_cv,function(x) max(1e-4,x)),qtrunc_cv=sapply(qtrunc_cv,function(x) max(1e-4,x)))%>%
  ggplot(aes(x=qmis_cv,y=qtrunc_cv,col=signif,label=gene_name))+
  geom_point()+
  scale_x_log10(breaks=10^(-4:0))+
  scale_y_log10(breaks=10^(-4:0))+
  ggrepel::geom_label_repel()+
  theme_classic()+
  geom_hline(yintercept=0.1,linetype=2)+
  geom_vline(xintercept=0.1,linetype=2)

#Combine individual dnds annotated muts output into single df
all_annotated_combined<-lapply(all_dndsout,function(dataset_dndsout) {
  data_set_annotmuts<-dplyr::bind_rows(lapply(dataset_dndsout,function(dndsout) dndsout$annotmuts))
  return(data_set_annotmuts)
  })%>%dplyr::bind_rows()

dndsout_combined_all_restricted$annotmuts%>%
  #filter(sampleID%in%published_samples)%>%
  filter(gene%in%all_ukbb_CH)%>%
  group_by(gene,impact)%>%
  summarise(n=n())%>%
  pivot_wider(names_from="impact",values_from = "n")%>%
  replace(is.na(.), 0)%>%
  print(n=30)

###TEST IF MUTATIONS ARE ON SHARED BRANCHES
ukbb_gene_muts<-dndsout_combined_published_restricted$annotmuts%>%
  filter(sampleID%in%published_samples)%>%
  filter(gene%in%novel_FI_ukbb & impact!="Synonymous")%>%
  dplyr::select(chr,pos,ref,mut)%>%
  apply(1,paste0,collapse="-")
length(ukbb_gene_muts)

pdf(paste("New_ukbb_shared_muts.pdf"),width=15,height=10)
Map(details=unlist(details_all_list,recursive = F),
    tree=unlist(trees_all_list,recursive=F),
    f=function(details,tree){
      details<-details%>%
        mutate(new_ukbb_shared=ifelse(mut_ref%in%ukbb_gene_muts & node>length(tree$tip.label),"yes","no"))
      
      if(any(details$new_ukbb_shared=="yes")){
        cat("New ukbb on shared branch found",sep="\n")
        #Plot the tree
        tree=plot_tree(tree,cex.label=0,vspace.reserve=0.05)
        temp=plot_tree_labels(tree,
                              details,
                              type="line",
                              query.field = "new_ukbb_shared", #alternative is 'coding_change_chip'
                              data.frame(value=c("yes"),col="red",pch = 17,stringsAsFactors = FALSE), #if use 'coding_change_chip', value is 'Coding change mutation in driver'
                              label.field = "Gene",
                              cex.label = 0.8,
                              lty=2,
                              lwd=2)
      }
    })
dev.off()

###Run dNdS by dataset
dnds_by_dataset_all_genes<-lapply(all_combined_by_dataset,function(dataset_details) {
  dndsout_dataset<-dndscv_on_details(details=dataset_details,
                                 outp=3)
  return(dndsout_dataset)
})

dnds_by_dataset_all_ukbb<-Map(dataset_details=all_combined_by_dataset,dataset=names(all_combined_by_dataset),f=function(dataset_details,dataset) {
  cat(dataset,sep="\n")
  if(dataset%in%c("MF","SDS")) {
    return(NULL)
  } else {
    dndsout_dataset<-dndscv_on_details(details=dataset_details,
                                       gene_list=all_ukbb_CH,
                                       outp=3)
    return(dndsout_dataset)
  }
})

dnds_by_dataset_novel_ukbb<-Map(dataset_details=all_combined_by_dataset,dataset=names(all_combined_by_dataset),f=function(dataset_details,dataset) {
  cat(dataset,sep="\n")
  if(dataset%in%c("MF","SDS")) {
    return(NULL)
  } else {
    dndsout_dataset<-dndscv_on_details(details=dataset_details,
                                       gene_list=novel_FI_ukbb,
                                       outp=3)
    return(dndsout_dataset)
  }
})

names(dnds_by_dataset_all_genes)=names(dnds_by_dataset_all_ukbb)=names(dnds_by_dataset_novel_ukbb)<-names(all_combined_by_dataset)

lapply(dnds_by_dataset_all_genes,function(dndscv_out) {
  dndscv_out$sel_cv%>%
    filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1)%>%
    #arrange(qglobal_cv)%>%
    pull(gene_name)
})

dnds_by_dataset_all_genes$SDS$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1)%>%
  #arrange(qglobal_cv)%>%
  pull(gene_name)


lapply(dnds_by_dataset_all_ukbb,function(dndscv_out) {
  if(!is.null(dndscv_out)) {
    dndscv_out$sel_cv%>%
      filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1)%>%
      #arrange(qglobal_cv)%>%
      pull(gene_name)
  } else {
      NULL
    }
})

lapply(dnds_by_dataset_novel_ukbb,function(dndscv_out) {
  if(!is.null(dndscv_out)) {
    dndscv_out$sel_cv%>%
      filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1)%>%
      #arrange(qglobal_cv)%>%
      pull(gene_name)
  } else {
    NULL
  }
})

dnds_by_dataset_all_genes$DK$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1)

##Emily's published data only
#Lymphoid
dndsout_Ly<-dndscv_on_details(details=all_combined_by_dataset$Ly,
                              gene_list=all_ukbb_CH,
                              outp=3)

#Lymphoid
dndsout_Ly<-dndscv_on_details(details=all_combined_by_dataset$Ly,
                               gene_list=all_ukbb_CH,
                               outp=3)
dndsout_Ly$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)

#CML
dndsout_CML<-dndscv_on_details(details=all_combined_by_dataset$TN,
                              #gene_list=novel_FI_ukbb,
                              outp=3)
dndsout_CML$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)

#General older adult haemopoiesis datasets
dndsout_myeloid<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[c("EM","NW","MSC_BMT","TN","CML","MF")]),
                               #gene_list=novel_FI_ukbb,
                               outp=3)
signif<-dndsout_myeloid$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)%>%
  pull(gene_name)

dndsout_myeloid$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)%>%
  
dndsout_myeloid_restricted<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[c("EM","NW","MSC_BMT","TN","CML","MF")])%>%
                                                filter(sampleID%in%published_samples),
                                   gene_list=novel_FI_ukbb,
                                   outp=3)
dndsout_myeloid_restricted$sel_cv%>%
  #filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)%>%
  head(n=20)

dndsout_myeloid$annotmuts%>%
  #filter(sampleID%in%published_samples)%>%
  filter(gene%in%signif)%>%
  group_by(sampleID,gene,impact)%>%
  summarise(n=n())%>%
  pivot_wider(names_from="impact",values_from = "n")%>%
  replace(is.na(.), 0)%>%
  print(n=30)

#NOW EXCLUDING THE UNPUBLISHED DATA
dndsout<-dndscv_on_details(details=combined_all%>%
                                     filter(sampleID%in%published_samples),
                                   outp=3)
signif<-dndsout_myeloid$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)%>%
  pull(gene_name)

dndsout$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)

dndsout_restricted<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[c("EM","NW","MSC_BMT","TN","CML","MF")])%>%
                                                filter(sampleID%in%published_samples),
                                              gene_list=novel_FI_ukbb,
                                              outp=3)
dndsout_restricted$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)

##########

dndsout_myeloid<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[c("EM","NW","MSC_BMT","TN","CML","MF")])%>%
                                     filter(sampleID%in%published_samples),
                                   #gene_list=novel_FI_ukbb,
                                   outp=3)
signif<-dndsout_myeloid$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)%>%
  pull(gene_name)

dndsout_myeloid$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)
  
dndsout_myeloid_restricted<-dndscv_on_details(details=dplyr::bind_rows(all_combined_by_dataset[c("EM","NW","MSC_BMT","TN","CML","MF")])%>%
                                                  filter(sampleID%in%published_samples),
                                                gene_list=all_ukbb_CH,
                                                outp=3)
dndsout_myeloid_restricted$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)

##########

dndsout_EMpublished<-dndscv_on_details(details=all_combined_by_dataset$NW%>%
                                     filter(sampleID%in%published_samples),
                                   #gene_list=novel_FI_ukbb,
                                   outp=3)

signif<-dndsout_EMpublished$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)%>%
  pull(gene_name)

dndsout_EMpublished$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)%>%slice_head(n=10)

dndsout_EMpublished_restricted<-dndscv_on_details(details=all_combined_by_dataset$EM%>%
                                                filter(sampleID%in%published_samples),
                                              gene_list=novel_FI_ukbb,
                                              outp=3)
dndsout_EMpublished_restricted$sel_cv%>%
  filter(qmis_cv<0.1|qtrunc_cv<0.1|qallsubs_cv<0.1|qglobal_cv<0.1)%>%
  arrange(qglobal_cv)


dndsout_myeloid$annotmuts%>%
  #filter(sampleID%in%published_samples)%>%
  filter(gene%in%signif)%>%
  group_by(sampleID,gene,impact)%>%
  summarise(n=n())%>%
  pivot_wider(names_from="impact",values_from = "n")%>%
  replace(is.na(.), 0)%>%
  print(n=30)