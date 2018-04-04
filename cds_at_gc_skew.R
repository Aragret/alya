rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

user = 'Alina'
if (user == 'Alina') {setwd('/home/aragret/Alina/other_projects/alya/')}

library(seqinr)
################################################################################
### genes.txt is a result of create_cds_table.py script, which takes a directory with 
### CDs and makes a table "species - gene - sequence"

df = read.table('genes.txt', sep = '\t')
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG')

he <- c()
for(i in 1:nrow(df)){
  # i = 1
  s <-  s2c(as.character(df[i, 3]))     #string into vector of chars
  Codons <- splitseq(s, 0, 3)   #vector of chars into codons
  CodonsTable = data.frame(Codons)
  CodonsTable2 = CodonsTable[CodonsTable$Codons %in% VecOfSynFourFoldDegenerateSites,]
  four_fold_number = length(CodonsTable2)
  CodonsTable2 = paste(CodonsTable2, collapse="")
  if(CodonsTable2 == ''){
    next
  }
  s <-  s2c(CodonsTable2) 
  Codons <- splitseq(s, 0, 3)
  nuctable = c()
  for (j in 1:length(Codons)){
    nuctable = rbind( nuctable, s2c(Codons [j])[3])
  } 
  counts <- as.data.frame(table(nuctable))
  str(counts)
  vec <- c("A", "T", "G", "C")
  counts <- counts[counts$nuctable %in% vec,]
  a = counts[counts$nuctable == 'A',]$Freq; t = counts[counts$nuctable == 'T',]$Freq
  g = counts[counts$nuctable == 'G',]$Freq; c = counts[counts$nuctable == 'C',]$Freq
  at_skew = (a - t) / (a + t)
  gc_skew = (g - c) / (g + c)
  he <- rbind(he, c(as.character(df[i, 1]), as.character(df[i,2]), at_skew, gc_skew, four_fold_number))
  
}
he = as.data.frame(he)
names(he) = c('Species', 'Gene', 'ATSkew', 'GCSkew', 'Sites_numbers')
write.table(he, file='cds_at_gc_skew_with_sites_numbers.txt', quote = FALSE, sep = '\t', row.names = FALSE)

summary(he[he$Gene == 'ATP8', 'ATSkew'])

hist(as.numeric(as.character(he[he$Gene == 'ATP8', 'ATSkew'])))
hist(as.numeric(as.character(he[he$Gene == 'ATP8', 'GCSkew'])))

hist(as.numeric(as.character(he[he$Gene == 'CytB', 'ATSkew'])))
hist(as.numeric(as.character(he[he$Gene == 'CytB', 'GCSkew'])))
