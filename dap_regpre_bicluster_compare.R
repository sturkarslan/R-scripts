# This script collects members of cMonkey biclusters, DAP-Chip targets and RegPrecise
# regulator targets, compare them each othe by using hypergeometric distribution
# and identifies common regulatory themes.
# by Serdar Turkarslan

######################## Variables ##########################
output.folder = "Desktop/Serdar_motifGetter/dap_chip_compare/"
genome = read.csv("Desktop/Serdar_motifGetter/dap_chip_compare/dvh_img_gene_info.txt", header=T, sep="\t")
annot = read.csv("Dropbox/dvh/DvH_Annotations.txt", header=T, sep="\t")
net.data =read.csv("~/Dropbox/dvh/DvHmotifNetwork.sif", sep="\t")
regprecise.short = read.csv("Desktop/Serdar_motifGetter/dap_chip_compare/regPrecise_short.txt", header=T, sep=",")
    colnames(net.data) = c("gene", "interaction", "bicluster")
dir = list.files(output.folder, pattern="fasta.meme")
inf.data =read.csv("Dropbox/dvh/infNetwork.sif", header=F, sep="\t")
    colnames(inf.data) = c("regulator", "interaction", "bicluster")
bicluster.d1 =read.csv("~/Dropbox/dvh/DvHmotifNetwork.sif", sep="\t")
  colnames(bicluster.d1) = c("gene", "interaction", "bicluster")
bicluster.enrichment = read.csv("Dropbox/dvh/bicluster_functional_enrichment.csv", header=T, sep=",")
regulon.properties = read.csv("Desktop/Serdar_motifGetter/dap_chip_compare/regPrecise_regulons.txt", header=T, sep=",")
dapchip.hits = data.frame(NULL)

#convert non DVU gene names for regprecise targets to DVU format
merged.0 = merge(regprecise.short, annot, by.x ="genes", by.y = "name", all.x=F, all.y=F) # long.list from regprecise_parsing.R
merged.1 = merged.0[,1:2]
merged.2 = merged.0[,10:10]
regprecise.data = cbind(merged.1, merged.2)
names(regprecise.data) = c("short.genes", "regulators", "genes")

# Combine inf regulator targets with bicluster member genes so 
# that you will have targets for each inf regulator
bc2inf = merge(bicluster.d1, inf.data, all.x=F, all.y=F, by="bicluster")

# select targets for each regulator in inferelator ( not considering AND gates)
# needed to cleanup motifs and AND gates
bc2inf.nomotif = bc2inf[grep("motif", bc2inf[,"gene"], invert=T),]
bc2inf.noand = bc2inf.nomotif[grep("AND", bc2inf.nomotif[,"regulator"], invert=T),]

# Collect dap-chip and ortholog motif sequences
for (i in dir) {
    filename = paste(output.folder, i, "/", substr(i, 1, 7), "_dap_orth.csv", sep="")
      if (file.exists(filename)) {
      cat("Processing: ", i, "\n")
      f1 = read.csv(filename, header=T, sep=",")
      f2 = cbind( Response.Reg = substr(i, start = 1, stop = 7),
                  Target.Id = as.character(f1[,"locusId"]),
                  Motif.Hit = f1[,"hit"],
                  Method = as.character(f1[,"method"]))
      dapchip.hits = rbind(dapchip.hits, f2)
       }
}
# Collect dapchip identified targets for  response regulator
dapchip.members = list()
    for(dapchip.i in unique(dapchip.hits[,"Response.Reg"])) {
      dapchip.members[[dapchip.i]] = dapchip.hits[which(dapchip.hits[,"Response.Reg"]==dapchip.i),"Target.Id"]
    }
# collect members of biclusters
bicluster.members = list()
    for(bicluster.i in unique(net.data[,"bicluster"])) {
      bicluster.members[[bicluster.i]] = net.data[which(net.data[,"bicluster"]==bicluster.i),"gene"]
    }
# Select RegPrecise regulon members for each regulator
 regprecise.targets = list()
    for(regprecise.i in unique(regprecise.data[,"regulators"])) {
      regprecise.targets[[regprecise.i]] = regprecise.data[which(regprecise.data[,"regulators"]==regprecise.i),"genes"]
    }
# Select Inferelator regulon members for each regulator
inf.gene.members = list()
for(inf.gene.i in unique(bc2inf.noand[,"regulator"])) {
  inf.gene.members[[inf.gene.i]] = bc2inf.noand[which(bc2inf.noand[,"regulator"]==inf.gene.i),"gene"]
}

# Create new background model for hypergeometric
Ncm.dc = intersect(unique(net.data[,"gene"]), unique(dapchip.hits[,"Target.Id"] ))
Ncm.rp = intersect(unique(net.data[,"gene"]), unique(regprecise.data[,"genes"] ))
Ndc.rp = intersect(unique(dapchip.hits[,"Target.Id"]), unique(regprecise.data[,"genes"] ))
Nif.rp = intersect(unique(bc2inf.noand[,"gene"]), unique(regprecise.data[,"genes"] ))

############# Enriched biclusters for dap.chip targets #############
pValues = matrix(ncol=length(dapchip.members),
                nrow=length(bicluster.members),
                dimnames= list(names(bicluster.members),
                names(dapchip.members)))
for(dapchip.j in names(dapchip.members)) {
  for(bicluster.j in names(bicluster.members)) {
    # paramaters for hypergeometric
    m.cm.dc = length(intersect(Ncm.dc, dapchip.members[[dapchip.j]]))
    k.cm.dc = length(intersect(Ncm.dc, bicluster.members[[bicluster.j]]))
    n.cm.dc = length(Ncm.dc) - m.cm.dc
    # hypergeometric p values
        pValues[bicluster.j, dapchip.j] = phyper(length(intersect(bicluster.members[[bicluster.j]], dapchip.members[[dapchip.j]])), # q
                                            m.cm.dc , # m
                                            n.cm.dc , # n
                                            k.cm.dc , # k
                                            lower.tail=F)
      
  }
}
  p.adj.Values = p.adjust(pValues, method = "BH")  # Benjamini & Hochberg p. values
  PP = (p.adj.Values > 0) & (p.adj.Values <= 1.0e-02)
  unique.dap.member = length(unique(dapchip.hits[,"Response.Reg"]))
  unique.bicluster.number = length(unique(net.data[,"bicluster"]))
  RN = matrix(rep(rownames(p.adj.Values), length(colnames(p.adj.Values))), nrow=unique.bicluster.number)
  CN = matrix(rep(colnames(p.adj.Values), length(rownames(p.adj.Values))), ncol=unique.dap.member, byrow=T)
  dapchip_biclus.BH = data.frame(bicluster=RN[PP], Response.Reg=CN[PP], BH.p.value=p.adj.Values[PP])
  write.csv(dapchip_biclus.BH, "Desktop/Serdar_motifGetter/dapchip_biclus_BH.csv")
######### Enriched regprecise regulons for dap.chip targets ##########
pValues = matrix(ncol=length(dapchip.members),
                nrow=length(regprecise.targets),
                dimnames= list(names(regprecise.targets),
                names(dapchip.members)))
for(dapchip.j in names(dapchip.members)) {
  for(regpre.j in names(regprecise.targets)) {
    # paramaters for hypergeometric
    m.dc.rp = length(intersect(Ndc.rp, dapchip.members[[dapchip.j]]))
    k.dc.rp = length(intersect(Ndc.rp, regprecise.targets[[regpre.j]]))
    n.dc.rp = length(Ndc.rp) - m.dc.rp
        # hypergeometric p values
        pValues[regpre.j, dapchip.j] = phyper(length(intersect(regprecise.targets[[regpre.j]], dapchip.members[[dapchip.j]])), # q
                                            m.dc.rp , # m
                                            n.dc.rp , # n
                                            k.dc.rp , # k
                                            lower.tail=F)
      
  }
}
  p.adj.Values = p.adjust(pValues, method = "BH")
  PP = (p.adj.Values > 0) & (p.adj.Values <= 1.0e-02)
  unique.dap.member = length(unique(dapchip.hits[,"Response.Reg"]))
  unique.regpre.number = length(unique(regprecise.data[,"regulators"]))
  RN = matrix(rep(rownames(p.adj.Values), length(colnames(p.adj.Values))), nrow=unique.regpre.number)
  CN = matrix(rep(colnames(p.adj.Values), length(rownames(p.adj.Values))), ncol=unique.dap.member, byrow=T)
  dapchip_regpre.BH = data.frame(bicluster=RN[PP], Response.Reg=CN[PP], BH.p.value=p.adj.Values[PP])

############# Enriched biclusters for regprecise targets #############
pValues = matrix(ncol=length(regprecise.targets),nrow=length(bicluster.members), dimnames= list(names(bicluster.members), names(regprecise.targets)))
for(reg.j in names(regprecise.targets)) {
  for(bicluster.j in names(bicluster.members)) {
    # paramaters for hypergeometric
    m.cm.rp = length(intersect(Ncm.rp, regprecise.targets[[reg.j]]))
    k.cm.rp = length(intersect(Ncm.rp, bicluster.members[[bicluster.j]]))
    n.cm.rp = length(Ncm.rp) - m.cm.rp
    # hypergeometric p values
        pValues[bicluster.j, reg.j] = phyper(length(intersect(bicluster.members[[bicluster.j]],regprecise.targets[[reg.j]])), # q
                                            m.cm.rp , # m
                                            n.cm.rp , # n
                                            k.cm.rp , # k
                                            lower.tail=F)
      
  }
}
  # Benjamini Hochberg correction  
  p.adj.Values = p.adjust(pValues, method = "BH") # Benjamini & Hochberg p. values
  PP = (p.adj.Values > 0) & (p.adj.Values <= 1.0e-04)
  unique.regpre.number = length(unique(regprecise.data[,"regulators"]))
  unique.bicluster.number = length(unique(net.data[,"bicluster"]))
  RN = matrix(rep(rownames(p.adj.Values), length(colnames(p.adj.Values))), nrow=unique.bicluster.number)
  CN = matrix(rep(colnames(p.adj.Values), length(rownames(p.adj.Values))), ncol=unique.regpre.number, byrow=T)
  regpre_biclus.BH = data.frame(bicluster=RN[PP], Regulator=CN[PP], BH.p.value=p.adj.Values[PP])
  merge.1 = merge(regpre_biclus.BH, bicluster.enrichment, by.x= "bicluster", by.y= "bicluster", all.x=T, all.y=F, sort=F)
  merge.2 = merge(merge.1, regulon.properties, by.x= "Regulator", by.y= "regulator", all.x=T, all.y=F, sort=F)
  write.csv(merge.2, "Desktop/Serdar_motifGetter/regpre_biclus_BH.csv")
  write.csv(regpre_biclus.BH, "Desktop/Serdar_motifGetter/regpre_biclus_BH.csv")
############# Enriched Inferelator targets for regprecise targets #############
pValues = matrix(ncol=length(regprecise.targets),nrow=length(inf.gene.members), dimnames= list(names(inf.gene.members), names(regprecise.targets)))
for(reg.j in names(regprecise.targets)) {
  for(inf.j in names(inf.gene.members)) {
    # paramaters for hypergeometric
    m.if.rp = length(intersect(Nif.rp, regprecise.targets[[reg.j]]))
    k.if.rp = length(intersect(Nif.rp, inf.gene.members[[inf.j]]))
    n.if.rp = length(Nif.rp) - m.if.rp
    # hypergeometric p values
        pValues[inf.j, reg.j] = phyper(length(intersect(inf.gene.members[[inf.j]],regprecise.targets[[reg.j]])), # q
                                            m.if.rp , # m
                                            n.if.rp , # n
                                            k.if.rp , # k
                                            lower.tail=F)
      
  }
}
  # Benjamini Hochberg correction  
  p.adj.Values = p.adjust(pValues, method = "BH") # Benjamini & Hochberg p. values
  PP = (p.adj.Values > 0) & (p.adj.Values <= 1.0e-04)
  unique.regpre.number = length(unique(regprecise.data[,"regulators"]))
  unique.inf.number = length(unique(bc2inf.noand[,"regulator"]))
  RN = matrix(rep(rownames(p.adj.Values), length(colnames(p.adj.Values))), nrow=unique.inf.number)
  CN = matrix(rep(colnames(p.adj.Values), length(rownames(p.adj.Values))), ncol=unique.regpre.number, byrow=T)
  regpre_inf.BH = data.frame(Inf.Regulator=RN[PP], RegPre.Regulator=CN[PP], BH.p.value=p.adj.Values[PP])
  write.csv(regpre_inf.BH, "Desktop/Serdar_motifGetter/regpre_inferelator_BH.csv")

# DVU0942 network
d1 = regprecise.data[which(regprecise.data[,"regulators"]=="DVU0942"),]
d2 = cbind(gene = as.character(d1[,"genes"]), interaction = c("regulates"), bicluster = as.character(d1[,"regulators"]))
e1 = net.data.1[which((net.data.1[,"bicluster"]=="bicluster_0034") | (net.data.1[,"bicluster"]=="bicluster_0155") | (net.data.1[,"bicluster"]=="bicluster_0342")| (net.data.1[,"bicluster"]=="bicluster_0151")),]
f1 = rbind(d2, e1)

# DVU0606 network
d3 = regprecise.data[which(regprecise.data[,"regulators"]=="DVU0606"),]
d4 = cbind(gene = as.character(d3[,"genes"]), interaction = c("regulates"), bicluster = as.character(d3[,"regulators"]))
e2 = net.data.1[which((net.data.1[,"bicluster"]=="bicluster_0206") | (net.data.1[,"bicluster"]=="bicluster_0081")),]
f2 = rbind(d4, e2)
write.csv(f2, "Dropbox/dvh/regprecise_0606.csv", quote=F)

# DVU3229 network
d3 = regprecise.data[which(regprecise.data[,"regulators"]=="DVU3229"),]
d4 = cbind(gene = as.character(d3[,"genes"]), interaction = c("regulates"), bicluster = as.character(d3[,"regulators"]))
e2 = net.data.1[which((net.data.1[,"bicluster"]=="bicluster_0122")),]
f2 = rbind(d4, e2)
write.csv(f2, "Dropbox/dvh/regprecise_3229.csv", quote=F)

# DVU0057 network
d1 = regprecise.data[which(regprecise.data[,"regulators"]=="DVU0057"),]
d2 = cbind(gene = as.character(d1[,"genes"]), interaction = c("regulates"), bicluster = as.character(d1[,"regulators"]))
e1 = net.data.1[which((net.data.1[,"bicluster"]=="bicluster_0005") | (net.data.1[,"bicluster"]=="bicluster_0316") | (net.data.1[,"bicluster"]=="bicluster_0094")),]
f1 = rbind(d2, e1)
write.csv(f1, "Dropbox/dvh/regprecise_0057.csv", quote=F)

## Venn Diagrams
layout(matrix(c(1,2,3), byrow = TRUE))
bs = c("bicluster_0034", "bicluster_0155", "bicluster_0342")
for (r in bs) {
  b1 = bicluster.members[[r]]
  r1 = regprecise.targets[["DVU0942"]]
  union = union(bicluster.members[[r]], regprecise.targets[["DVU0942"]])
  Counts <- matrix(0, nrow=length(union), ncol=2)
  colnames(Counts) <- c(names(bicluster.members[r]), names(regprecise.targets["DVU0942"]))
for (i in 1:length(union))
  {
  Counts[i,1] <- union[i] %in% b1
  Counts[i,2] <- union[i] %in% r1
  }
vennDiagram( vennCounts(Counts), counts.col = c("blue"), circle.col = c("orange", "green"))

  }

bs = c("bicluster_0034", "bicluster_0155", "bicluster_0342")

  b1 = bicluster.members[["bicluster_0034"]]
  b2 = bicluster.members[["bicluster_0155"]]
  b3 = bicluster.members[["bicluster_0342"]]
  r1 = regprecise.targets[["DVU0942"]]
  union = union(regprecise.targets[["DVU0942"]], union(bicluster.members[["bicluster_0034"]], bicluster.members[["bicluster_0155"]]))
  union = union[grep("motif", union, invert=T)]
  Counts <- matrix(0, nrow=length(union), ncol=3)
  colnames(Counts) <- c(names(bicluster.members["bicluster_0034"]),  names(bicluster.members["bicluster_0155"]), names(regprecise.targets["DVU0942"]))
for (i in 1:length(union))
  {
  Counts[i,1] <- union[i] %in% b1
  Counts[i,2] <- union[i] %in% b2
  Counts[i,3] <- union[i] %in% r1
  }
vennDiagram( vennCounts(Counts), counts.col = c("blue"), circle.col = c("orange", "green",  "red"))

