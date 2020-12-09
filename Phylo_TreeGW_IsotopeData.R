library(V.PhyloMaker)
library(Taxonstand)
library(diversitree)
library(ape)
library(caper)
library(phytools)

#Read in Fan et al (2017) plant rooting database
infile = "E:/Global_Rooting/Plant_Iso_List.csv"
RD = read.csv(infile,header=TRUE)
RD$FullName = paste(RD$Genus,RD$Species)

#Validate plant names against "The Plant List"
r1 <- TPL(RD$FullName, corr = TRUE)
TPL_Names = data.frame(species = r1$Taxon, genus = r1$Genus, family = r1$Family)
PlantNames = TPL_Names[which(TPL_Names$family != ""),]
PlantNames$matchname = gsub(" ", "_", PlantNames$species)

# Run Phyolomaker
result <- phylo.maker(PlantNames, scenarios=c("S1","S2","S3"))
result_analysis = result$scenario.3
result_analysis$node.label <- NULL


#Assign Iso as traits
for (i in 1:length(result_analysis$tip.label))
{
  speciesname = result_analysis$tip.label[i]
  family = as.character(PlantNames[PlantNames$matchname == speciesname,3])
  Records = RD[paste(RD$Genus,"_",RD$Species,sep="") == speciesname,]
  maxIso = round(max(Records$Iso))    
  result_analysis$trait.state_Iso[i] = maxIso
  result_analysis$matchname[i] = paste(Records$Genus[1],"_",Records$Species[1],sep="")
  result_analysis$family[i] = family[1]
}

result_export = result_analysis
result_export$tip.label = paste(result_export$family,".",result_export$matchname,sep="")
write.tree(result_export, file = "E:/Global_Rooting/Phylo_Newick/ISO_Newick.tree")

Export_Annotation = data.frame(result_analysis$family, result_analysis$tip.label, result_analysis$trait.state_Iso)
write.table(Export_Annotation,"E:/Global_Rooting/Phylo_Newick/ISO_Newick_Annotation.tree",sep=" ",row.names = FALSE,col.names = FALSE)

#Reformat traits
Traits = data.frame(result_analysis$trait.state_Iso)
row.names(Traits) <- result_analysis$tip.label
colnames(Traits) <- c("Iso")
write.csv(Traits,"Iso_Evaristo.csv")

#D-test for trait occurrence
Traits_Dtest = data.frame(result_analysis$matchname, result_analysis$trait.state_Iso)
colnames(Traits_Dtest) <- c("matchname","Iso")
Trees = comparative.data(result_analysis, Traits_Dtest, matchname)
redPhyloD <- phylo.d(Trees, binvar=Iso)
print(redPhyloD)

#Phylogeny plot
tiff("E:/Global_Rooting/Figures/Phylogeny_ISO_031820.tiff", units="in", width=5, height=5, res=600)
trait.plot(result_analysis, dat=Traits, class=result_analysis$family, cols = list(Iso = c("lightblue1","blue4")),cex.lab=0.4)
text(x=0, y=60, paste("n = ",length(result_analysis$tip.label)," species",sep=""),col ="black", cex=0.8)
text(x=0, y=30, paste("D = ",round(redPhyloD$DEstimate,3), " (p =",round(redPhyloD$Pval1,3),")",sep=""),col ="blue", cex=0.8)
dev.off()

