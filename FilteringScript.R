library(tidyverse)

#UV Presence Absence 
UV <- read.delim("/Users/paulwaters/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_3_Rsx/Submission_GenomeBiology/R1/Manuscript_revised/Resubmission/proteinGroups_UV.txt") %>% 
  filter((Peptides.UV.test.1 > 1 | Peptides.UV.test.2 > 1) &
           Peptides.UV.control.1 == 0 &
           Peptides.UV.control.2 == 0 &
           Peptides.UV.control.3 == 0) %>% 
  select(Protein.IDs, Peptides.UV.test.1, Peptides.UV.test.2, Peptides.UV.control.1, Peptides.UV.control.2, Peptides.UV.control.3,Potential.contaminant,Reverse) %>% 
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+")

#PFA Presence Absence
PFA <- read.delim("/Users/paulwaters/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_3_Rsx/Submission_GenomeBiology/R1/Manuscript_revised/Resubmission/proteinGroups_PFA.txt") %>% 
  filter((Peptides.ALL > 1 | Peptides.PFA.all > 1) &
  Peptides.PFA.none == 0) %>%
  select(Protein.IDs, Peptides.ALL, Peptides.PFA.all, Peptides.O3, Peptides.PFA.none, Potential.contaminant,Reverse) %>% 
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+")

#Native Presence Absence
Native <- read.delim("/Users/paulwaters/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_3_Rsx/Submission_GenomeBiology/R1/Manuscript_revised/Resubmission/proteinGroups_native.txt") %>% 
  filter(rowSums(cbind(Peptides.2 > 1, Peptides.4 > 1, Peptides.5 > 1, Peptides.6 > 1, Peptides.7 > 1, Peptides.ALL > 1, Peptides.combine > 1)) >= 2 &
           Peptides.Scram == 0 &
           Peptides.N1 == 0 &
           Peptides.N2 == 0 &
           Peptides.none == 0)  %>% 
  select(Protein.IDs, Peptides.2, Peptides.4, Peptides.5, Peptides.6, Peptides.7, Peptides.ALL, Peptides.combine, Peptides.Scram, Peptides.N1, Peptides.N2, Peptides.none, Potential.contaminant,Reverse) %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+")

All_Pres_Abs <- bind_rows(UV %>% select(1), PFA %>% select(1), Native %>% select(1)) %>%
  distinct()


#UV Ratio
UV_R <- read.delim("/Users/paulwaters/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_3_Rsx/Submission_GenomeBiology/R1/Manuscript_revised/Resubmission/proteinGroups_UV.txt") %>% 
  filter(Peptides.UV.control.1 != 0 | Peptides.UV.control.2 != 0 | Peptides.UV.control.3 != 0) %>% 
  select(Protein.IDs, LFQ.intensity.UV.control.1, LFQ.intensity.UV.control.2, LFQ.intensity.UV.test.1, LFQ.intensity.UV.test.2, Potential.contaminant,Reverse) %>% 
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>% 
  mutate(UV_Ratio = log2((LFQ.intensity.UV.test.1 + LFQ.intensity.UV.test.2) / (LFQ.intensity.UV.control.1 + LFQ.intensity.UV.control.2))) %>% 
  filter(UV_Ratio > 1.584963, UV_Ratio != Inf)

#PFA Ratio
PFA_R <- read.delim("/Users/paulwaters/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_3_Rsx/Submission_GenomeBiology/R1/Manuscript_revised/Resubmission/proteinGroups_PFA.txt") %>% 
  filter(Peptides.PFA.none != 0) %>% 
  select(Protein.IDs, LFQ.intensity.PFA.none, LFQ.intensity.PFA.all, LFQ.intensity.O3, LFQ.intensity.ALL, Potential.contaminant,Reverse) %>% 
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>% 
  mutate(PFA_Ratio = log2((((LFQ.intensity.PFA.all + LFQ.intensity.ALL)/2) / LFQ.intensity.PFA.none))) %>% 
  mutate(PFA_O3 = log2((LFQ.intensity.O3 / LFQ.intensity.PFA.none))) %>%
  filter(PFA_Ratio != Inf, PFA_O3 != Inf & (PFA_Ratio > 1.584963 | PFA_O3 > 1.584963))

#Native Ratio
Native_R <- read.delim("/Users/paulwaters/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_3_Rsx/Submission_GenomeBiology/R1/Manuscript_revised/Resubmission/proteinGroups_native.txt") %>%
  mutate(across(starts_with("LFQ.intensity"), as.numeric)) %>%
  mutate(Native_Ratio = log2(((LFQ.intensity.2 + LFQ.intensity.4 + LFQ.intensity.5 + LFQ.intensity.6 + LFQ.intensity.7 + LFQ.intensity.ALL + LFQ.intensity.O3)/7) / ((LFQ.intensity.N1 + LFQ.intensity.N2 + LFQ.intensity.none)/3))) %>% 
  filter(Native_Ratio > 1.584963) %>% 
  filter(rowSums(cbind(LFQ.intensity.2 > 1, LFQ.intensity.4 > 1, LFQ.intensity.5 > 1, LFQ.intensity.6 > 1, LFQ.intensity.7 > 1, LFQ.intensity.ALL > 1, LFQ.intensity.O3 > 1)) > 4) %>%
  select(Protein.IDs, LFQ.intensity.N1, LFQ.intensity.N2, LFQ.intensity.none, LFQ.intensity.2, LFQ.intensity.4, LFQ.intensity.5, LFQ.intensity.6, LFQ.intensity.7, LFQ.intensity.ALL, LFQ.intensity.O3, Potential.contaminant, Reverse)
  #filter(Potential.contaminant != "+" & Reverse != "+") %>%

Native_R$pvalue_all <- sapply(1:nrow(Native_R), function(i) t.test(Native_R[i,c("LFQ.intensity.2","LFQ.intensity.4","LFQ.intensity.5","LFQ.intensity.6","LFQ.intensity.7","LFQ.intensity.ALL")], Native_R[i,c("LFQ.intensity.N1", "LFQ.intensity.N2","LFQ.intensity.none")],var.equal = T,alternative = "greater")$p.value)

Native_R$p.adj <- p.adjust(Native_R$pvalue_all, method = "BH")

Native_R <- Native_R %>% 
  filter(p.adj < 0.05)

All_Ratios <- bind_rows(UV_R %>% select(1), PFA_R %>% select(1), Native_R %>% select(1)) %>% 
  distinct()

all <- bind_rows(All_Pres_Abs, All_Ratios) %>% 
  distinct()

all <- all %>%
  separate(col = Protein.IDs, into = c("Protein.ID1", "Protein.ID2", "Protein.ID3"), sep = "\\|") %>%
  distinct(Protein.ID1) 

write.table(all, file = "/Users/paulwaters/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_3_Rsx/Submission_GenomeBiology/R1/Manuscript_revised/Resubmission/All_Pres_Abs_Ratios.txt", sep = "\t", quote = F, row.names = F)

