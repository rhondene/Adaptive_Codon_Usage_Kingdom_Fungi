## Phylogenetic analysis for Figure 1 data
## Author: Rhondene Wint
#PIC between GC3 and Codon usage
library(ape)
library(phangorn)
library(seqinr)
library(phytools)
library(ggtree)
library(treeio)
library(geiger)
library(nlme)
#library(dplyr)
library(caper)



fungal_tree = read.tree("./all_codon_fracsUpdated.newick")
is.rooted(fungal_tree)
is.binary(fungal_tree) #binary and rooted -->


#load codon bias and individual codon freq trait data : I.e. codfreq/amino acid usage

cod_freq = read.csv('./all_codon_frac_traitsUpdated.csv', row.names="Species") 
GC3 = cod_freq[,'Mean_GC3']
ENC = cod_freq[,'Mean_ENC']
names(GC3) <- names(ENC) <- rownames(cod_freq)

#### PIC for ENC~GC3
GC3_pic = pic(GC3, fungal_tree)
ENC_pic = pic(ENC, fungal_tree)
pic_model=lm(ENC_pic~GC3_pic - 1)
summary(pic_model)


####### GC3,ENC Maximum likelihood modefitting
trait = ENC
bm = fitContinuous(phy=fungal_tree, dat=trait, model='BM')
ou = fitContinuous(phy=fungal_tree, dat=trait, model='OU')
eb = fitContinuous(phy=fungal_tree, dat=trait, model='EB')
del = fitContinuous(phy=fungal_tree, dat=trait, model='delta')

ou$opt$aicc
bm$opt$aicc
eb$opt$aicc
del$opt$aicc

#gc3
trait = GC3
bm = fitContinuous(phy=fungal_tree, dat=trait, model='BM')
ou = fitContinuous(phy=fungal_tree, dat=trait, model='OU')
eb = fitContinuous(phy=fungal_tree, dat=trait, model='EB')
del = fitContinuous(phy=fungal_tree, dat=trait, model='delta')

ou$opt$aicc
bm$opt$aicc
eb$opt$aicc
del$opt$aicc

####---------compute PIC for each codon 

#load tatbe rscu just to ge tthe codon names
rscu_traits = read.csv('./all_rscu_traits.csv', row.names="Species") 
codons = names(rscu_traits)
# store each PIC results
pic_summary=list()

for (codon in codons){
  cod=cod_freq[,codon]
  cod_pic = pic(cod,fungal_tree)
  pic_model=lm(GC3_pic~cod_pic - 1)
  codon_pic=summary(pic_model)
  #add to list
  pic_summary[[codon]]=codon_pic
}



pic_summary2=list()
p_vals= list()

for (codon in codons){
  cod=cod_freq[,codon]
  cod_pic = pic(cod,fungal_tree)
  pic_model=lm(GC3_pic~cod_pic - 1)
  codon_pic=summary(pic_model)
  #add to list
  pic_summary2[[codon]]=c(codon_pic$coefficients,codon_pic$adj.r.squared)
}

write.csv(pic_summary2,"./codfreq_gc3_v2.csv")


####-----do PIC (enc~codfreq)
cod_freq = read.csv('./all_codon_frac_traits.csv', row.names="Species") 
enc = cod_freq[,'Mean_ENC']
ENC_pic = pic(enc, fungal_tree)

pic_summary_enc=list()

for (codon in codons){
  cod=cod_freq[,codon]
  cod_pic = pic(cod,fungal_tree)
  pic_model=lm(ENC_pic~cod_pic - 1)
  codon_pic=summary(pic_model)
  #add to list
  pic_summary_enc[[codon]]=c(codon_pic$coefficients,codon_pic$adj.r.squared)
}

write.csv(pic_summary_enc,"./codfreq_enc.csv")

##do PHYLOGENETIC CORRECTED Pearson correlation between each codon vs ENC and GC3,separately
####----------------------------- do pic(enc~codfreq)

pic_summary_enc=list()

for (codon in codons){
  cod=cod_freq[,codon]
  cod_pic = pic(cod,fungal_tree)
  res=cor.test(ENC_pic,cod_pic, method='pearson')
  #add to list
  pic_summary_enc[[codon]]=c(res$estimate,res$p.value)
}

write.csv(pic_summary_enc,"./PIC_pearson_codfreq_enc.csv")

##do GC3
pic_summary_gc=list()

for (codon in codons){
  cod=cod_freq[,codon]
  cod_pic = pic(cod,fungal_tree)
  res=cor.test(GC3_pic,cod_pic, method='pearson')
  #add to list
  pic_summary_gc[[codon]]=c(res$estimate,res$p.value)
}

write.csv(pic_summary_gc,"./PIC_pearson_codfreq_gc.csv")

