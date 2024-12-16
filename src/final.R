if(interactive()) setwd('~/src/vpspr_collab/') 
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(openxlsx))
suppressMessages(library(smoother))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(plotrix))
suppressMessages(library(LDlinkR))
suppressMessages(library(corrplot))
suppressMessages(library(gtools))
suppressMessages(library(grid))
suppressMessages(library(arrow))
suppressMessages(library(binom))
select    = dplyr::select
summarize = dplyr::summarize


# Function ####
fp = function(x1, n1, x2, n2) {
  fisher.test(matrix(c(x1, n1-x1, x2, n2-x2), nrow=2, byrow=T), simulate.p.value=T)$p.value
}
fe = function(x1, n1, x2, n2) {
  fisher.test(matrix(c(x1, n1-x1, x2, n2-x2), nrow=2, byrow=T), simulate.p.value=T)$estimate
}


find_max_ac = function(af,an,ci=.95) {
  if (af == 0) {
    return (0)
  } else {
    quantile_limit = ci # ci for one-sided, 1-(1-ci)/2 for two-sided
    max_ac = qpois(quantile_limit,an*af)
    return (max_ac)
  }
}

find_af_filter = function(ac, an, ci=.95, lower=(.1/(max_an)), upper=2, tol=1e-7, precision=1e-6) { 
  if (is.na(ac) | is.na(an) | ac == 0 | an == 0 | ac == 1) {
    return (0.0)
  } else {
    quantile_limit = ci 
    attempt_uniroot = tryCatch({
      uniroot_result = uniroot(f = function(af,ac,an) { return (ac - 1 - qpois(p=quantile_limit,lambda=an*af)) },lower=lower,upper=upper,ac=ac,an=an,tol=tol)
    }, warning = function(w) {
      print(paste("ac= ",as.character(ac),", an= ",as.character(an)," warning = ",as.character(w),sep=''))
      return (0.0)
    }, error = function(e) {
      print(paste("ac= ",as.character(ac),", an= ",as.character(an)," error = ",as.character(e),sep=''))
      return (0.0)
    }, finally = {
    })
    max_af = round(uniroot_result$root,-log10(precision)) # round to nearest millionth
    while(find_max_ac(af=max_af,an=an,ci=ci) < ac) {
      max_af = max_af + precision # increase by millionths until you find the upper bound - the highest AF for which 95%CI AC is still less than observed AC
    }
    max_af = max_af - precision # back off one unit from the AF that violated the 95%CI AC < obs AC condition
    return (max_af)
  }
}

find_af_filter_vectorized = function(ac_vector, an_vector, ci=.95, lower=(.1/(max_an)), upper=2, tol=1e-7) { 
  return (mapply(find_af_filter, ac_vector, an_vector, ci=ci, lower=lower, upper=upper, tol=tol))
}

meansd = function(x, digits=1) {
  paste0(formatC(mean(x, na.rm=T),format='f',digits=digits),'±',formatC(sd(x, na.rm=T),format='f',digits=digits))
}

tla_to_ola = function(x) {
  mapping = data.frame(tla=c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","TER","THR","TRP","TYR","VAL"),
                       ola=c("A",  "R",  "N",  "D",  "C",  "Q",  "E",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",  "S",  "X",  "T",  "W",  "Y",  "V"))
  for (row in 1:dim(mapping)[1]) {
    x = gsub(mapping$tla[row], mapping$ola[row], toupper(x))
  }
  return (x)
}

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

# creat table 2####
demo = read_tsv("data/analytic/vpspr_samples.tsv") %>%
  mutate(age_at_death = suppressWarnings(as.integer(age_at_death)))

demo %>%
  group_by(gt129) %>%
  summarize(.groups='keep', 
            n = n(),
            mean_age = mean(age_at_death, na.rm=T),
            sd_age   =   sd(age_at_death, na.rm=T),
            meansd_age = meansd(age_at_death,digits=0),
            mean_dur = mean(duration_months, na.rm=T),
            sd_dur   =   sd(duration_months, na.rm=T),
            meansd_dur = meansd(duration_months/12),
            sex = paste0(sum(sex=='M'),'M/',sum(sex=='F'),'F'),
            n_pos_hx = sum(family_hx_neurological_disease %in% 'Yes'),
            n_neg_hx = sum(family_hx_neurological_disease %in% 'No'),
            p_hx = sum(family_hx_neurological_disease %in% 'Yes')/sum(family_hx_neurological_disease %in% c('Yes','No'))) %>%
  ungroup() %>%
  mutate(hx = paste0(percent(p_hx,digits=0),' (',n_pos_hx,'/',n_neg_hx+n_pos_hx,')')) -> demo_smry


demo %>%
  summarize(.groups='keep', 
            n = n(),
            mean_age = mean(age_at_death, na.rm=T),
            sd_age   =   sd(age_at_death, na.rm=T),
            meansd_age = meansd(age_at_death,digits=0),
            mean_dur = mean(duration_months, na.rm=T),
            sd_dur   =   sd(duration_months, na.rm=T),
            meansd_dur = meansd(duration_months/12),
            sex = paste0(sum(sex=='M'),'M/',sum(sex=='F'),'F'),
            n_pos_hx = sum(family_hx_neurological_disease %in% 'Yes'),
            n_neg_hx = sum(family_hx_neurological_disease %in% 'No'),
            p_hx = sum(family_hx_neurological_disease %in% 'Yes')/sum(family_hx_neurological_disease %in% c('Yes','No'))) %>%
  ungroup() %>%
  mutate(hx = paste0(percent(p_hx,digits=0),' (',n_pos_hx,'/',n_neg_hx+n_pos_hx,')')) %>%
  mutate(gt129='all') %>% 
  relocate(gt129) -> demo_total_row



demo_smry_disp = rbind(demo_smry, demo_total_row) %>%
  select("Codon 129 genotype"= gt129, 
         N=n, 
        "Age at death(years,mean±SD)"= meansd_age, 
        "Disease duration (years,mean±SD)" = meansd_dur, 
        Sex=sex, 
        "Family history ofneurological disease"=hx)


write_tsv(demo_smry_disp,"display_items/table-2.tsv")



age_anova = aov(age_at_death ~ gt129, data=demo)
summary(age_anova)
TukeyHSD(age_anova)

dur_anova = aov(duration_months ~ gt129, data=demo)
summary(dur_anova)
TukeyHSD(dur_anova)

binom.test(x=sum(demo$sex=="M"), n=sum(!is.na(demo$sex)),p=0.5)
sex_age_anova = aov(age_at_death ~ sex, data=demo)
sex_dur_anova =aov(duration_months ~ sex, data=demo)
summary(sex_age_anova)
summary(sex_dur_anova)



sex_fisher = fisher.test(table(demo[,c('sex','gt129')]))



output_file <- "output/statistics_results.txt"

sink(output_file)

cat("AGE ANOVA RESULTS:\n")
print(summary(age_anova))
cat("\nPost-hoc Tukey Test for AGE ANOVA:\n")
print(TukeyHSD(age_anova))

cat("\nDURATION ANOVA RESULTS:\n")
print(summary(dur_anova))
cat("\nPost-hoc Tukey Test for DURATION ANOVA:\n")
print(TukeyHSD(dur_anova))

cat("\nSEX AGE ANOVA RESULTS:\n")
print(summary(sex_age_anova))

cat("\nSEX DURATION ANOVA RESULTS:\n")
summary(sex_dur_anova)


cat("\nFISHER'S EXACT TEST RESULTS:\n")
print(sex_fisher)

sink()

cat("Results saved to:", output_file, "\n")






# MERGE ALL DATASET####

## constants ####

m129v_pos = 4699605
tibial_nerve_eqtl_lead_pos = 4663722
cerebellum_eqtl_lead_pos = 4613886
e200k_haplo_marker_pos = 4611775
# position 4611775 (rs1046265426) is a variant found exclusively on E200K haplotypes in our gPrD dataset
# (all 29/29 individuals het for this variant are also het for E200K)

filtered_color = 'gray'
linked_color = 'purple'
antilinked_color = 'green'
m129v_color = 'orange'
linked_r2_threshold = 0.25
default_color = 'black'
eqtl_color = '#46a4ad'
vpspr_col= "#00C5CD"


original_threshold = 4.3e-6
max_credible_af = original_threshold
loose_threshold = original_threshold*5

## sequencing depth #### 
gene_depth = read_tsv("data/analytic/gene_depth_smry.tsv", col_types=cols()) 


## gnomAD DATA ####
gnomad_prnp = read.table('data/analytic/gnomAD_v4.1.0_20-4570000-4709000_2024_05_24_10_57_32.csv', sep = ",", header = TRUE) %>%
  as_tibble() %>%
  clean_names()  %>%
  dplyr::rename(chrom = chromosome, 
                pos = position,
                ref = reference, 
                alt = alternate) %>%
  dplyr::mutate(chrom = as.character(chrom)) %>%
  mutate(pos_id = paste(chrom, formatC(pos, width=9, flag='0', format='f', digits=0), ref, alt, sep='-'))

## SYNTHETIC VCF - annotation of PRNP variants
# synthetic VCF with some annotations
#synth = read.table('reference/prnp_synthetic.table',sep='\t',header=T)
#colnames(synth) = gsub('[^a-z0-9_]','_',tolower(colnames(synth)))
grch38_coding_offset = 4699605 - 4680251 # 19354. based on rs1799990
synth = read_tsv('data/analytic/prnp_synthetic.tsv', col_types=cols()) %>%
  clean_names() %>%
  mutate(pos = pos + grch38_coding_offset) %>%
  mutate(pos_id = paste(chrom, formatC(pos, width=9, flag='0', format='f', digits=0), ref, alt, sep='-')) %>%
  mutate(codon = (cds_position - 1) %/% 3 + 1) %>%
  mutate(syn = nchar(amino_acids)==1 & !is.na(amino_acids))

# get annotations easier to read
hgvsp_split = strsplit(synth$hgv_sp,',')
synth$hgvsp1 = mapply('[[',hgvsp_split,1)
hgvsp1_change = strsplit(synth$hgvsp1,'\\.')
synth$hgvsp1_change = mapply('[',hgvsp1_change,3)
synth$hgvsp1_change[grepl('>',synth$hgvsp1_change)] = NA
synth$aa_change = tla_to_ola(synth$hgvsp1_change)




##TARGET SUMMARY####

targeted_wide =read_csv("data/analytic/targeted_summary.csv", col_types=cols())

targeted_annotated = targeted_wide %>%
  left_join(gnomad_prnp %>% select(pos_id, group_max_faf_group, group_max_faf_frequency,ref,alt), by = "pos_id")


position_list = unique(targeted_annotated$pos_id)





## ODDS RATIO####
###gprd ####

pseud0 = 0.05
pseudinf_gprd = 20

odds_ratio_gprd <- targeted_annotated %>%
  filter(
    is.finite(ac_VPSPr), !is.na(ac_VPSPr),
    is.finite(an_VPSPr), !is.na(an_VPSPr),
    is.finite(ac_gPrD), !is.na(ac_gPrD),
    is.finite(an_gPrD), !is.na(an_gPrD)
  ) %>%
  rowwise() %>%  
  mutate(fisher_p_gprd = fp(ac_VPSPr, an_VPSPr, ac_gPrD, an_gPrD)) %>%
  mutate(or_gprd = case_when(ac_VPSPr==0 & ac_gPrD==0 ~ 1, # replace NaN with 1 when both are absent
                                TRUE ~ fe(ac_VPSPr, an_VPSPr, ac_gPrD, an_gPrD))) %>%
  mutate(or_gprd_pseudo=pmin(pmax(or_gprd,pseud0),pseudinf_gprd)) %>%
  ungroup() %>%
  select(pos_id,or_gprd,fisher_p_gprd,or_gprd_pseudo) %>%
  right_join(targeted_annotated, by ="pos_id" ) 

odds_ratio_gprd$p_bonf_gprd = pmin(1, odds_ratio_gprd$fisher_p_gprd * nrow(odds_ratio_gprd))



###NFE ####

NFE_gnomAD <- gnomad_prnp %>%
  inner_join(targeted_wide, by='pos_id') %>%
  select(pos_id,ref, alt,ac_nfe= allele_count_european_non_finnish,an_nfe= allele_number_european_non_finnish,ac_VPSPr,an_VPSPr,filter)

pseud0 = 0.05
pseudinf_nfe = 20

odds_ratio_nfe <- NFE_gnomAD %>%
  filter(
    is.finite(ac_VPSPr), !is.na(ac_VPSPr),
    is.finite(an_VPSPr), !is.na(an_VPSPr),
    is.finite(ac_nfe), !is.na(ac_nfe),
    is.finite(an_nfe), !is.na(an_nfe)
  ) %>%
  rowwise() %>%  
  mutate(fisher_p_nfe = fp(ac_VPSPr, an_VPSPr, ac_nfe, an_nfe)) %>%
  mutate(or_nfe = case_when(ac_VPSPr==0 & ac_nfe==0 ~ 1, # replace NaN with 1 when both are absent
                                    TRUE ~ fe(ac_VPSPr, an_VPSPr, ac_nfe, an_nfe))) %>%
  mutate(or_nfe_pseudo=pmin(pmax(or_nfe,pseud0),pseudinf_nfe)) %>%
  ungroup() %>%
  select(pos_id,or_nfe,fisher_p_nfe,or_nfe_pseudo) %>%
  right_join(NFE_gnomAD, by ="pos_id" )

odds_ratio_nfe$p_bonf_nfe = pmin(1, odds_ratio_nfe$fisher_p_nfe * nrow(odds_ratio_nfe))




###sCJD_GRCh37####
sCJD_GRCh37= read_tsv("data/analytic/sCJD_GRCh37.tsv", col_types=cols())

sCJD_GRCh38= read_tsv("data/analytic/sCJD_GRCh37_liftover_to_GRCh38.tsv", col_types=cols()) %>%
  select(base_pair_location,position=base_pair_location_hg38)%>%
  inner_join(sCJD_GRCh37,by='base_pair_location') %>%
  select(-base_pair_location)%>%
  mutate(pos_id=paste0(chromosome, "-00", position, "-", other_allele, "-", effect_allele))%>%
  select(last_col(),everything())


read_csv("data/analytic/filter_pos.csv",col_names =FALSE, show_col_types = FALSE) %>% ### just to match to r^2 exported from hail, r^2 value are exported without pos_id
  t() %>%
  as_tibble(rownames='index') %>%
  mutate(index = as.integer(gsub('X','',index))) %>%
  mutate(pos_id =  gsub("'","",V1)) %>% # remove single quotes
  mutate(pos_id = gsub('20-','20-00',pos_id)) %>% # convert to 9-digit pos 
  select(index, pos_id) -> filter_pos



###  r_squared ####
r_squared = read_tsv("data/analytic/r_squared.tsv",col_names =FALSE, show_col_types = FALSE) 

colnames(r_squared) = filter_pos$pos_id
r_squared <- r_squared %>%
  mutate(pos_id = filter_pos$pos_id) %>%
  select(last_col(), everything())

r_squared <- as.data.frame(r_squared)
row.names(r_squared)=r_squared$pos_id

r_squared %>%
  select(-pos_id) -> r_squared

r_squared[upper.tri(r_squared)] <- t(r_squared)[upper.tri(r_squared)]


r_squared <- r_squared %>%
  mutate(pos_id = row.names(r_squared)) %>%
  select(last_col(), everything())

row.names(r_squared) <- 1:nrow(r_squared)



variant_4699605 = r_squared %>%
  filter(pos_id == "20-004699605-A-G") %>%
  select(-pos_id) %>%
  t() %>%
  as_tibble() %>%
  mutate(r_squared = as.numeric(V1)) %>%
  mutate(pos_id = filter_pos$pos_id) %>%
  select(pos_id, r_squared)

variant_4692886 = r_squared %>%
  filter(pos_id == "20-004692886-G-A") %>%
  select(-pos_id) %>%
  t() %>%
  as_tibble() %>%
  mutate(r_squared = as.numeric(V1)) %>%
  mutate(pos_id = filter_pos$pos_id) %>%
  select(pos_id, r_squared)




odds_ratio_gprd %>%
  left_join(variant_4699605, by='pos_id') %>%
  relocate(pos_id_orig) %>%
  mutate(group_max_faf_frequency = replace_na(group_max_faf_frequency, 0)) %>%
  mutate(group_max_faf_group = case_when(is.na(group_max_faf_group) | group_max_faf_group=='' ~ 'none',
                                         TRUE ~ group_max_faf_group)) -> targeted_all




max_an = max(targeted_all$an_gPrD)
targeted_all$faf_gprd = find_af_filter_vectorized(targeted_all$ac_gPrD, targeted_all$an_gPrD)


pseud0= 1e-06
targeted_all <- targeted_all %>%
  mutate(faf_combined = pmax(group_max_faf_frequency,pseud0,faf_gprd))


## PARQUET FROM GTEx ####

### nerve_tibial ####
prnp_tibial = read_tsv("data/analytic/gtex_prnp_tibial_eqtl.tsv", col_types=cols())



### cerebellum ####
prnp_cerebellum = read_tsv("data/analytic/gtex_prnp_cerebellum_eqtl.tsv", col_types=cols())







all_dataset= targeted_all %>%
  select(pos_id_orig,pos_id,qual,filter,faf_gprd,group_max_faf_group,group_max_faf_frequency,faf_combined)%>%
  left_join(odds_ratio_gprd %>% 
              select(pos_id,
                     or_gprd,
                     fisher_p_gprd,
                     ac_VPSPr,an_VPSPr,ac_gPrD,an_gPrD,
                     p_bonf_gprd,
                     or_gprd_pseudo),by='pos_id')%>%
  left_join(odds_ratio_nfe %>%
              select(pos_id,
                     or_nfe,
                     fisher_p_nfe,
                     ac_nfe,an_nfe,
                     p_bonf_nfe,
                     or_nfe_pseudo),by='pos_id') %>%
  left_join(sCJD_GRCh38 %>%
              select(pos_id,
                     af_scjd_gwas=effect_allele_frequency,
                     or_scjd_gwas=odds_ratio,
                     pval_scjd_gwas=p_value),by='pos_id') %>%
  left_join(variant_4699605 %>%
              select(pos_id,
                     r_squared_4699605=r_squared),by='pos_id')%>%
  left_join(variant_4692886 %>%
              select(pos_id,
                     r_squared_4692886=r_squared),by='pos_id')%>%
  left_join(prnp_tibial%>%
              select(pos_id,
                     slope_tibial_nerve=slope,
                     pval_nominal_tibial_nerve=pval_nominal),by='pos_id')%>%
  left_join(prnp_cerebellum%>%
              select(pos_id,
                     slope_cerebellum=slope,
                     pval_nominal_cerebellum=pval_nominal),by='pos_id')%>%
  mutate(postion= str_extract(pos_id, "(?<=-)[0-9]+(?=-)"),
         postion = str_replace(postion, "^0+", ""),
         pos_id_1=pos_id) %>%
  separate(pos_id_1, into = c(NA, NA, "ref", "alt"), sep = "-")%>%
  select(pos_id_orig,pos_id,postion,ref,alt,everything())









##EXOME DATA ####

faf95 = read_tsv('data/analytic/exome_sites.tsv.gz', col_types=cols())


pseud0=1e-6

## plot pmax
faf95 <- faf95 %>%
  mutate(af_pseudo = pmax(faf95_gnomad,incon_faf95,pseud0))


#### get gene length

constraint_metrics = read_tsv("data/analytic/gnomad.v4.1.constraint_metrics.tsv", col_types=cols())

max_cds_length =  constraint_metrics %>%
  select(gene, cds_length) %>%
  filter(!is.na(cds_length)) %>%
  arrange(desc(cds_length)) %>%
  filter(!duplicated(gene)) %>%
  rename(symbol=gene)


vep_unique = read_tsv("data/analytic/exome_vep.tsv") 



# filter for impact, AF, and presence in a gene
vep_unique_filtered = vep_unique %>%
  filter(symbol != '-') %>%
  filter(impact %in% c('MODERATE','HIGH')) %>%
  left_join(faf95, by="pos_id") %>%
  filter(faf95_gnomad <= original_threshold & incon_faf95 <= original_threshold) 

# group by gene
vep_gby_gene_orig_step1 = vep_unique_filtered %>%
  group_by(symbol) %>%
  summarise(sum_vpspr_ac=sum(vpspr_ac),
            .groups = "keep") %>%
  ungroup()

vep_gby_gene_orig_threshold = max_cds_length %>%
  inner_join(vep_gby_gene_orig_step1, by="symbol")




vep_gby_gene_loose_step1 = vep_unique %>%
  filter(symbol != '-') %>%
  filter(impact %in% c('MODERATE','HIGH')) %>%
  left_join(faf95, by="pos_id") %>%
  filter(faf95_gnomad <= loose_threshold & incon_faf95 <= loose_threshold) %>%  
  group_by(symbol) %>%
  summarise(sum_vpspr_ac=sum(vpspr_ac),
            .groups = "keep") %>%
  ungroup()

vep_gby_gene_loose_threshold = max_cds_length %>%
  inner_join(vep_gby_gene_loose_step1, by="symbol")




#PLOT PANEL####



## Figure 1 ####

prevalance_all_prion_disease= 1/6239 # Maddox 2020
vpspr_proportion_prion_disease = 88/3931 # https://case.edu/medicine/pathology/research/national-prion-disease-pathology-surveillance-center/surveillance/tables-cases-examined Accessed May 24, 2024, last updated April 23, 2024; there is a May 19 snapshot on web.archive.org
penetrance = 0.42 # proportion with positive family history from Notari 2018,  15/36 confirmed by Brian
one_in = 1/(prevalance_all_prion_disease*vpspr_proportion_prion_disease)
max_credible_af = 4.3e-6 # from calculator, https://www.cardiodb.org/allelefrequencyapp/ see Whiffin & Minikel 2017

orhx=read_tsv("data/analytic/country_case_smry.tsv",show_col_types = FALSE)
gnomad_prnp_plot=read_tsv("data/analytic/population_control_af.tsv",show_col_types = FALSE)


resx=300
png('display_items/figure-1.png',width=6.5*resx,height=3.5*resx,res=resx)
par(mfrow=c(1,2))
panel = 1

par(mar=c(3,4,3,1))
xlims = c(0, 1)
ylims = c(1, 100000)
yats = rep(1:9, 5) * 10^rep(0:4, each=9)
ybigs = 10^(0:4)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F, log='y')
axis(side=1, at=0:20/20, labels=NA, tck=-0.025)
axis(side=1, at=0:10/10, labels=NA, tck=-0.05)
axis(side=1, at=0:10/10, lwd=0, line=-0.5, labels=percent(0:10/10), cex.axis=0.8)
mtext(side=1, line=1.5, text='proportion of cases with family history', cex=0.8)
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, labels=formatC(ybigs,format='d',big.mark=','), lwd=0, line=-0.5, las=2, cex.axis=0.8)
mtext(side=2, line=2.75, text='odds ratio', cex=0.8)
points(x=orhx$hx, y=orhx$or, pch=19)
segments(x0=orhx$hx, y0=orhx$or_l95, y1=orhx$or_u95, lwd=1.5)
segments(x0=orhx$hx_l95, x1=orhx$hx_u95, y0=orhx$or, lwd=1.5)
par(xpd=T)
#text(x=orhx$hx, y=orhx$or, labels=orhx$aa_change, srt=45, pos=4)
text(x=orhx$hx, y=orhx$or*.5, labels=orhx$aa_change, pos=4, cex=0.7)
par(xpd=F)
x = 0:100/100
m = loess(log10(or) ~ hx, data=orhx, span=1.5)
y = predict(m, x)
points(x, 10^y, type='l', lty=3)
famhx = binom.confint(x=15, n=36, method='wilson') # 15/36 in Notari 2018, Brian confirmed
abline(v=famhx$mean, col=vpspr_col, lty=1)
mtext(side=3, line=0, at=famhx$mean, col=vpspr_col, text='VPSPr', cex=0.8)
rect(xleft=famhx$lower, xright=famhx$upper, ybottom=min(ylims), ytop=max(ylims), border=NA, col=alpha(vpspr_col,0.1))
abline(v=famhx[,c('lower','upper')], col=vpspr_col, lty=3)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1



#vpspr_or_point_estimate = 10^predict(m, famhx_notari2018$mean)


par(mar=c(3,4,3,1))
ac_case = c(1:67) # range of possible VPSPr ACs assuming dominant inheritance
n_case = 3931 # total U.S. surveillance prion cases since 2008, regardless of subtype
n_control = max(gnomad_prnp_plot$n_control) # 807077

xlims = c(1/(2*n_control), 1e-2)
ylims = c(0, 67)
xbigs = c(10^(-7:-1), .5)
xbiglabs = c('1e-7','1e-6', '1e-5', '1e-4', '0.1%', '1%', '10%', '50%')
xats = rep(1:9, 7) * rep(10^(-7:-1), each=9)

plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i',log='x')
axis(side=1, at=xbigs, tck=-0.03, labels=NA)
axis(side=1, at=xbigs, lwd=0, labels=xbiglabs, line=-0.5, cex.axis=0.8)
axis(side=1, at=xats, labels=NA, tck=-0.015)
mtext(side=1, line=1.6, text='population control allele frequency',cex=0.8)
axis(side=2, at=0:5*30,las=2, cex.axis=0.8)
axis(side=2, at=0:10*10, labels=NA,tck=-0.02)
mtext(side=2, line=2.25, text='VPSPr AC', cex=0.8)

abline(v=max_credible_af, col='red', lty=3)
mtext(side=3, line=0.25, at=max_credible_af, col='red', text='max\ncredible\nAF', cex=0.8)

or = 1000
ac_control = ac_case * n_control / (or * n_case)
af_control = ac_control / n_control
points(x=ac_control/n_control, y=ac_case, type='l', lty=3, lwd=0.5)
mtext(side=3, at=max(af_control[ac_case < max(ylims)]), line=0.25, text='OR=1K', cex=0.7)

polygon(x=c(af_control[af_control < max_credible_af], max_credible_af, max_credible_af, rev(af_control[af_control < max_credible_af])),
        y = c(ac_case[af_control < max_credible_af], max(ac_case[af_control < max_credible_af]), 67, rev(rep(67, sum(af_control < max_credible_af)))),
        border=NA, col=alpha(vpspr_col, 0.1))

text(x=1.5e-6, y=40, col=vpspr_col, labels='hypoth-\nesized\nresults', cex=0.65)

or = 100
ac_control = ac_case * n_control / (or * n_case)
af_control = ac_control / n_control
points(x=ac_control/n_control, y=ac_case, type='l', lty=3, lwd=0.5)
mtext(side=3, at=max(af_control[ac_case < max(ylims)]), line=0.25, text='OR=100', cex=0.7)

or = 10
ac_control = ac_case * n_control / (or * n_case)
af_control = ac_control / n_control
points(x=ac_control/n_control, y=ac_case, type='l', lty=3, lwd=0.5)
mtext(side=3, at=max(af_control[ac_case < max(ylims)]), line=0.25, text='OR=10', cex=0.7)

mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1


dev.off()

## Figure 2: exomes  ####
if('figure-2'=='figure-2') {
  resx =300
  
  png('display_items/figure-2.png', width=6.5*resx, height=3*resx, res=resx)
  
  
  layout_matrix = matrix(c(1,1,1,2,2,2), nrow=1, byrow=T)
  
  layout(layout_matrix,
         widths=c(7,7,7,7,7,7),
         heights=14)
  
  #layout.show(7)
  
  
  pseud0=1e-6
  ac_case = c(1:134)
  or = 1000
  n_control =  730947 + 76215 # gnomad v4.1 counts
  n_case = 3931 # total U.S. surveillance prion cases since 2008, regardless of subtype
  ac_control = ac_case * n_control / (or * n_case)
  af_control = ac_control / n_control
  
  
  ###  exome af  ####
  
  par(mar=c(4,5,4,2))
  xlims = c(pseud0/5,1)
  ylims = c(1,135)
  xbigs = c(10^(-6:-1), .5)
  xbiglabs = c('0', '1e-5', '1e-4', '0.1%', '1%', '10%', '50%')
  xats = rep(1:9, 6) * rep(10^(-6:-1), each=9)
  xats = xats[3:length(xats)]
  
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F,log="x")
  axis(side=1, at=xlims, lwd.ticks=1,tck=-0.02,labels = NA)
  axis(side=1, at=xbigs, labels=NA)
  axis(side=1, at=xbigs, labels=xbiglabs, line=-0.5, lwd=0)
  axis.break(axis=1, breakpos=2e-6, style='slash')
  axis(side=1, at=xats, labels=NA, tck=-0.015)
  mtext(side=1, line=2, text='filtering allele frequency')
  axis(side=2, at=0:5*30, labels=NA, tck=-0.03)
  axis(side=2, at=0:5*30, las=2, line=-0.25, lwd=0)
  axis(side=2, at=0:10*10, labels=NA,tck=-0.015)
  mtext(side=2, line=2.2, text='VPSPr AC')
  subs = faf95 %>% filter(filter=='PASS')
  
  
  # Plot indels with gray 
  points(subs$af_pseudo[nchar(subs$ref) != nchar(subs$alt)], 
         subs$vpspr_ac[nchar(subs$ref) != nchar(subs$alt)], 
         cex = 0.1, col = "gray")
  
  # Plot SNPs black 
  points(subs$af_pseudo[nchar(subs$ref) == nchar(subs$alt)], 
         subs$vpspr_ac[nchar(subs$ref) == nchar(subs$alt)], 
         cex = 0.1, col = "black")

  
  legend(1e-5,130,  c('SNPs','Indels'),pch=20,col =c( 'black','gray'),cex=.8)
  
  #points(subs$af_pseudo, subs$vpspr_ac,cex=0.1)
  
  #faf95<- faf95 %>%
  #  filter(nchar(ref)==nchar(alt)) # only plot SNP
  
  subs = faf95 %>%
    filter(number_chrom=='chr20' & pos==m129v_pos)
  points(subs$af_pseudo, subs$vpspr_ac,pch=1,lwd=2,col=m129v_color)
  text(subs$af_pseudo, subs$vpspr_ac,col=m129v_color, pos=2, cex=0.8, labels=expression(paste(italic('PRNP'),' M129V')))
  polygon(x=c(af_control[af_control < max_credible_af], max_credible_af, max_credible_af, rev(af_control[af_control < max_credible_af])),
          y = c(ac_case[af_control < max_credible_af], max(ac_case[af_control < max_credible_af]), 134, rev(rep(134, sum(af_control < max_credible_af)))),
          border=NA, col=alpha(vpspr_col, 0.1))
  text(x=pseud0,y=120,col=vpspr_col,labels='expected\nhits',cex=0.8)
  abline(v=4.3e-6,col="red", lwd=1, lty=3 )
  mtext(text="4.3e-6", side=3, at=4.3e-6, line=0.25, col="red", cex=0.8)
  
  mtext("A", side=3, cex=1.5, adj = 0, line = 0.5)
  
  faf95_out = faf95 %>% 
    rename(chrom=number_chrom) %>% 
    mutate(pos_id = paste0(chrom,'-',formatC(pos,flag='0',width=9,format='d'),'-',ref,'-',alt)) %>%
    arrange(chrom, pos, ref, alt)
  write_tsv(faf95_out, 'output/exome_summary_stats.tsv.gz')
  
  # faf95 %>% 
  #   filter(filter=='PASS') %>% 
  #   filter(faf95_gnomad > 0.1 & vpspr_ac/132 > 2 * af_pseudo) %>% View()
  
  #rect(ybottom = 9, ytop = 132/2, xleft=3e-7, xright=3e-6, border='#00FFFF', col = NA, lwd=2)
  
  ### gene length ####
  
  par(mar=c(4,4,4,3))
  xlims = c(100,120000)
  ylims = c(0,67)
  xbigs = c(10^(2:5))
  xbiglabs = c('100bp','1kb','10kb','100kb')
  xats = rep(1:9, 6) * rep(10^(2:4), each=9)
  
  orig_color = '#000000'
  loose_color = '#C5E3BF'
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F,log="x")
  axis(side=1, at=xlims, lwd.ticks=1,tck=-0.02,labels = NA)
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, labels=xbiglabs, lwd=0, line=-0.25)
  axis(side=1, at=xats, labels=NA, tck=-0.02)
  mtext(side=1, line=2, text='CDS length')
  axis(side=2, at=0:5*30, labels=NA, tck=-0.03)
  axis(side=2, at=0:5*30, las=2, line=-0.25, lwd=0)
  axis(side=2, at=0:10*10, labels=NA,tck=-0.02)
  mtext(side=2, line=2.2, text='VPSPr AC')
  points(vep_gby_gene_loose_threshold$cds_length,vep_gby_gene_loose_threshold$sum_vpspr_ac,col=loose_color,pch=20,cex=0.8)  
  points(vep_gby_gene_orig_threshold$cds_length,vep_gby_gene_orig_threshold$sum_vpspr_ac,col=orig_color,pch=20,cex=0.8)
  
  par(xpd=T)
  legend(x=40000, y= max(ylims), c(original_threshold, loose_threshold),col =c(orig_color, loose_color),pch = 20,cex=0.8,title = "FAF threshold")
  mtext(text="AC=16", side=3, at=1e+05, line=-11, col="red", cex=0.8)
  par(xpd=F)
  
  abline(h = 16, lty = 2,col="red") 
  mtext("B", side=3, cex=1.5, adj = 0, line = 0.5)
  
  silence_is_golden = dev.off()
  
  vep_gby_gene_orig_threshold %>%
    left_join(vep_gby_gene_loose_threshold, by=c('symbol','cds_length'), suffix=c('_strict','_loose'))  -> sumbygene_out
  write_tsv(sumbygene_out, 'output/exome_summary_by_gene.tsv.gz')
}


## Figure 3 - targeted sequencing #### 
if('figure-3'=='figure-3') {
resx =300
png('display_items/figure-3.png', width=8*resx, height=5*resx, res=resx)




layout_matrix = matrix(c(1,1,1,3,3,3,
                         2,2,2,3,3,3,
                         4,4,4,5,5,5
), nrow=3, byrow=T)

layout(layout_matrix,
       widths=c(7,7,7,7,7,7),
       heights=c(10,4,14))




panel =1


## prnp position
human_fa = read.table('data/analytic/human_prnp_hg37_20_4663808_4682964.fa',skip=1)
human_sequence = paste(human_fa$V1, collapse='')
nchar(human_sequence)

human_fa_start = 4686456
human_fa_end = 4701590

exon1_start = c(4686456,4721909 )# transcription start site
exon1_end = c(4686512,4721969)

exon2_start = 4689138
exon2_end = 4689236


exon3_start =c(4699211,4724541)
exon3_end = c(4701588,4728460) # transcription end site

cds_start = c(4699221,4724552)
cds_end = c(4699982,4725079)




landmarks = tibble(gene = c('prnp','prnd'),
                   exon1_start = exon1_start,
                   exon2_start = exon2_start,
                   exon3_start = exon3_start,
                   exon1_end   = exon1_end,
                   exon2_end   =  exon2_end,
                   exon3_end   = exon3_end,
                   cds_start   = cds_start,
                   cds_end     =   cds_end,
                   exon2_fill =  '#A3A3A3',
                   y = 2)

### gene depth ####


xlims = c(4500000,4750000)
ylims = c(1,1250)
xbigs = seq(min(xlims),max(xlims),5e4)
xats = seq(min(xlims),max(xlims),1e4)

par(mar=c(1,6,2.5,3))

plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F,log="xy")
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, line=-0.5, lwd=0, labels=paste0(formatC(xbigs/1e6,format='f',digits=2), 'M'))
axis(side=1, at=xats, tck=-0.02, labels=NA)



par(las=2)
yats = rep(1:9, 4) * rep(10^(0:3), each=9)
yats = yats[3:length(yats)]
axis(side=2, at=yats, labels=NA, tck=-0.015)
ybigs = c(10^(0:3))
axis(side=2, at=ybigs, labels=ybigs)
mtext(side=2, line=2, text='Depth')
par(las=0)




polygon(x=c(gene_depth$posrounded,rev(gene_depth$posrounded)), y=c(gene_depth$l95,rev(gene_depth$u95)), col=alpha('#0000FF',0),border = alpha('#0000FF',0.2),lwd=2)
points(gene_depth$posrounded, gene_depth$avg_depth,type ="l",lwd=0.8)

abline(h=30,col="red", lwd=1, lty=3 )



mtext("A", side=3, cex=1.5, adj = 0, line = 0.5)

intron_lwd = 1
utr_lwd = 10
cds_lwd = 20

default_fill = '#000000'
utr_height = 2
cds_height = 4
prnd_height = 1





### annote prnp region ####

par(mar=c(2,6,2.5,3))

ylims = c(0, 2.5)
xlims = c(4500000,4750000)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
mtext(side=1, text='Position',line = 1.5,cex=1)


par(xpd=T)
segments(x0=landmarks$exon1_start, x1=landmarks$exon3_end, y0=landmarks$y, lwd=intron_lwd)
rect(xleft=landmarks$exon1_start, xright=landmarks$exon1_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=default_fill, border=NA)
#rect(xleft=landmarks$exon2_start, xright=landmarks$exon2_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=landmarks$exon2_fill, border='#000000')
rect(xleft=landmarks$exon3_start, xright=landmarks$exon3_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=default_fill, border=NA)
rect(xleft=landmarks$cds_start, xright=landmarks$cds_end, ybottom=landmarks$y-cds_height/2, ytop=landmarks$y+cds_height/2, col=default_fill, border=NA)

text(x=4694022, y=landmarks$y-3, labels="PRNP",cex=.7)

# PRND gene structure
text(x=4725184, y=landmarks$y-3, labels="PRND",cex=.7)


# annotate rs17327121
points(cerebellum_eqtl_lead_pos, 4, pch = 17, bg = "black")

text(cerebellum_eqtl_lead_pos, 4, pos=1, labels="rs17327121\n cerebellum\n eQTL lead SNP",cex=.6)


points(tibial_nerve_eqtl_lead_pos, 4, pch = 17, bg = "black")
text(tibial_nerve_eqtl_lead_pos, 4, pos=1, labels="rs6052749\n tibial nerve\n eQTL lead SNP",cex=.6)



scale_x = c(4750000,4755000)

arrows(x0=scale_x[1], x1=scale_x[2], y0=2, code=3, angle=90, length=0.005)
text(x=mean(scale_x), y=0.5, labels='5 kb',cex=0.7)



par(xpd=F)



###  targeted af ####
pseud0=1e-6


ac_case = c(1:134) 
n_case = 3931 # total U.S. surveillance prion cases since 2008, regardless of subtype
n_control =  730947 + 76215
max_credible_af= 4.3e-06

or = 1000
ac_control = ac_case * n_control / (or * n_case)
af_control = ac_control / n_control


par(mar=c(4,4,2.5,2))
xlims = c(pseud0/5,1)
ylims = c(0,135)

plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F,log="x")
axis(side=1, at=xlims, lwd.ticks=1,tck=-0.02,labels = NA)

xbigs = c(10^(-6:0))
xbiglabs = c('0','1e-5', '1e-4', '0.1%', '1%', '10%', '100%')
axis(side=1, at=xbigs, labels=xbiglabs)
axis.break(axis=1, breakpos=2e-6, style='slash')

xats = rep(1:9, 6) * rep(10^(-6:-1), each=9)
xats = xats[3:length(xats)]
axis(side=1, at=xats, labels=NA, tck=-0.015)
mtext(side=1, line=2.0, text='filtering allele frequency')

par(las=2)
axis(side=2, at=0:5*30)
axis(side=2, at=0:10*15, labels=NA,tck=-0.02)
par(las=0)
mtext(side=2, line=2.5, text='VPSPr AC')


points(pmax(all_dataset$faf_combined,pseud0), targeted_all$ac_VPSPr, pch=20, cex=0.3,
       col=case_when(all_dataset$filter != 'PASS' ~ filtered_color,
                     all_dataset$r_squared_4699605 >= linked_r2_threshold ~ linked_color, 
                     TRUE ~ default_color))

points(pmax(all_dataset$faf_combined,pseud0), targeted_all$ac_VPSPr, pch=20, cex=0.3,
       col=case_when(all_dataset$filter != 'PASS' ~ filtered_color,
                     all_dataset$r_squared_4699605 >= linked_r2_threshold ~ linked_color, 
                     TRUE ~ default_color))
points(all_dataset$faf_combined[all_dataset$postion==4699605],
       all_dataset$ac_VPSPr[all_dataset$postion==4699605],
       pch=1,col=m129v_color,lwd=0.5)
text(all_dataset$faf_combined[all_dataset$postion==4699605],
       all_dataset$ac_VPSPr[all_dataset$postion==4699605],
       pos=2, labels='129V',col=m129v_color,cex=0.8)

polygon(x=c(af_control[af_control < max_credible_af], max_credible_af, max_credible_af, rev(af_control[af_control < max_credible_af])),
        y = c(ac_case[af_control < max_credible_af], max(ac_case[af_control < max_credible_af]), 135, rev(rep(135, sum(af_control < max_credible_af)))),
        border=NA, col=alpha(vpspr_col, 0.1))

text(x=pseud0,y=120,col=vpspr_col,labels='expected\nhits',cex=0.8)
legend(x=1e-5, y=130, legend=c('PASS variants','filtered variants','r^2 ≥ 0.25 to 129V'),
       col=c(default_color, filtered_color, linked_color), pch=20, cex=0.8, bty='n')

abline(v=4.3e-6,col="red", lwd=1, lty=3 )
mtext(text="4.3e-6", side=3, at=4.3e-6, line=1, col="red", cex=0.8)

#xpected = 10^seq(-6,0,by=0.01)
#ypected = xpected * max(all_dataset$an_VPSPr)
#points(xpected, ypected, type='l', lwd=0.25)

mtext("B", side=3, cex=1.5, adj = 0, line = 0.5)





### gPrD vs. VPSPr ####


pseud0=1e-4

all_dataset_plot_af <- all_dataset %>%
  mutate(plot_af_gPrD = pmax(ac_gPrD/an_gPrD,pseud0),
         plot_af_VPSPr = pmax(ac_VPSPr/an_VPSPr,pseud0)
  )


par(mar=c(5,6,2,3))
xlims = c(pseud0/2,1)
ylims = c(pseud0/2,1)

plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F,log="xy")
axis(side=1, at=xlims, lwd.ticks=1,tck=-0.02,labels = NA)

xbigs = c(10^(-4:0))
xbiglabs = c('0' ,'0.1%', '1%', '10%', '100%')
axis(side=1, at=xbigs, labels=xbiglabs,cex.axis=0.8)
axis.break(axis=1, breakpos=2e-4, style='slash')

xats = rep(1:9, 4) * rep(10^(-4:-1), each=9)
xats = xats[3:length(xats)]
axis(side=1, at=xats, labels=NA, tck=-0.015)
mtext(side=1, line=2.5, text='gPrD AF')


axis(side=2, at=ylims, lwd.ticks=1,tck=-0.02,labels = NA)
par(las=2)
ybigs = c(10^(-4:0))
ybiglabs = c('0' ,'0.1%', '1%', '10%', '100%')
axis(side=2, at=ybigs, labels=ybiglabs, cex.axis=0.8)
axis.break(axis=2, breakpos=2e-4, style='slash')

yats = rep(1:9, 4) * rep(10^(-4:-1), each=9)
yats = yats[3:length(yats)]
axis(side=2, at=yats, labels=NA, tck=-0.015)

par(las=0)
mtext(side=2, line=3.5, text='VPSPr AF')


points(all_dataset_plot_af$plot_af_gPrD, all_dataset_plot_af$plot_af_VPSPr, pch=20, cex=0.5,
       col  = case_when(!is.na(all_dataset_plot_af$r_squared_4699605) & all_dataset_plot_af$r_squared_4699605 >= 0.25 ~ "purple",
                        !is.na(all_dataset_plot_af$r_squared_4699605) & all_dataset_plot_af$r_squared_4699605 <= -0.25 ~ "green",
                        TRUE ~ "grey"))
subs = all_dataset_plot_af %>% filter(postion == m129v_pos)
points(subs$plot_af_gPrD, subs$plot_af_VPSPr, pch=1, col=m129v_color, lwd=2)
par(xpd=T)
text(subs$plot_af_gPrD, subs$plot_af_VPSPr, col=m129v_color, pos=3, labels='129V')
par(xpd=F)

par(xpd=T)
legend(x=min(xlims)*1.5, y= max(ylims), c("r^2 ≥ 0.25",
                                          "-0.25 < r^2 < 0.25",
                                          "r^2 ≤ -0.25"),col =c( linked_color,'grey',antilinked_color),pch = 20,cex=0.8,bty='n')
par(xpd=F)

ors = c(0.01, 0.1, 10, 100)

line_color = "#00C5CD"

# correct
for (or in ors) {
  an_x = max(all_dataset_plot_af$an_gPrD)
  an_y = max(all_dataset_plot_af$an_VPSPr)
  ac_x = c(0.002 * an_x, 1:an_x)
  ref_x = an_x - ac_x
  ac_y = (or * an_y)/(or + ref_x / ac_x)
  ref_y = an_y - ac_y
  # only plot points with at least 1 ref or alt allele
  points_to_plot = ac_x > 1 & ac_x < an_x - 1 & ac_y > 1 & ac_y < an_y -1
  points(x=(ac_x/an_x)[points_to_plot], y=(ac_y/an_y)[points_to_plot], type='l', col=line_color, lty=3)
  label_x = min((ac_x/an_x)[points_to_plot])
  label_y = (ac_y/an_y)[which(ac_x/an_x >= label_x)[1]]
  text(label_x, label_y, srt=30, adj=0, col=line_color, labels=paste0('OR=',or), pos=ifelse(or > 1, 3, 1),cex=0.8)
}



mtext("C", side=3, cex=1.5, adj = 0, line = 0.5)





###  odds_ratio_gPrD #####


all_dataset_plot_or = all_dataset %>%
  filter(ac_VPSPr > 1 & ac_VPSPr < an_VPSPr - 1 & ac_gPrD > 1 & ac_gPrD < an_gPrD - 1 & ac_nfe/an_nfe > 0.001 & ac_nfe/an_nfe < 0.999) %>%
  filter(filter=='PASS') %>%
  filter(nchar(ref)==1 & nchar(alt)==1) %>%
  mutate(abs_or_gprd = case_when(or_gprd_pseudo > 1 ~ or_gprd_pseudo,
                                 or_gprd_pseudo < 1 ~ 1/or_gprd_pseudo)) %>%
  mutate(sign_or = or_gprd_pseudo > 1) %>%
  mutate(abs_or_nfe = case_when(or_nfe_pseudo > 1 ~ or_nfe_pseudo,
                                or_nfe_pseudo < 1 ~ 1/or_nfe_pseudo)) %>%
  mutate(conservative_or = case_when(sign_or ~ pmin(abs_or_gprd, abs_or_nfe, na.rm=T),
                                     !sign_or ~ pmax(1/abs_or_gprd, 1/abs_or_nfe, na.rm=T)))

M129V = all_dataset_plot_or[all_dataset_plot_or$pos_id == "20-004699605-A-G", ]
unpass = all_dataset_plot_or[all_dataset_plot_or$filter !="PASS", ]


leg <- data.frame(
  variant = c('M129V'),
  color = c('red')
)


pseud0=0.05

par(mar=c(5,4,2,2))
xlims = range(-1,1.15)
ylims = c(pseud0,pseudinf_gprd+10)

plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F,log="y")
axis(side=1, at=-4:2*0.5, lwd.ticks=1)
axis(side=1, at=-10:4*0.25,labels=NA,tck= -0.015)
mtext(side=1, line=2.5, text='r^2 to 129V')


par(las=2)
axis(side=2, at=c(pseud0/2,pseudinf_gprd+10), lwd.ticks=0,labels = NA)
ybigs = c(0.05,0.1,1,10, 20)
ybiglabs = c('0','0.1' ,'1', '10', '20')


axis(side=2, at=ybigs, labels=ybiglabs)
axis.break(axis=2, breakpos=0.06, style='slash')
par(las=0)

yats = rep(1:9, 4) * rep(10^(-2:1), each=9)
yats = yats[7:(length(yats)-7)]
axis(side=2, at=yats, labels=NA, tck=-0.015)
mtext(side=2, line=2.5, text='OR')




points(all_dataset_plot_or$r_squared_4699605, all_dataset_plot_or$conservative_or, pch=19,cex = 0.3) #case_when(all_dataset_plot_or$p_bonf_gprd < 0.01 ~ 1.5, all_dataset_plot_or$p_bonf_gprd >= 0.01 ~ 0.15))
points(M129V$r_squared_4699605, M129V$or_nfe, col=m129v_color, pch=1, lwd=2)
text(x=M129V$r_squared_4699605, y=M129V$or_nfe,labels='129V',pos=3,cex=0.8, col=m129v_color)
abline(v=0,col="red", lwd=1, lty=3 )
abline(h=1,col="red", lwd=1, lty=3 )


y_label_offset = 1.6
x_label_offset = 0.04
subs = all_dataset_plot_or %>% filter(postion==tibial_nerve_eqtl_lead_pos)
points(subs$r_squared_4699605, subs$conservative_or, pch=1, col=eqtl_color, lwd=2) #case_when(all_dataset_plot_or$p_bonf_gprd < 0.01 ~ 1.5, all_dataset_plot_or$p_bonf_gprd >= 0.01 ~ 0.15))
segments(x0=subs$r_squared_4699605, x1=subs$r_squared_4699605+x_label_offset, y0=subs$conservative_or, y1=subs$conservative_or/y_label_offset, col=eqtl_color)
text(subs$r_squared_4699605, subs$conservative_or/y_label_offset, col=eqtl_color, pos=4, labels='rs6052749 (tibial nerve eQTL lead SNP)', cex=0.8) 



subs = all_dataset_plot_or %>% filter(postion==cerebellum_eqtl_lead_pos)
points(subs$r_squared_4699605, subs$conservative_or, pch=1, col=eqtl_color, lwd=2) #case_when(all_dataset_plot_or$p_bonf_gprd < 0.01 ~ 1.5, all_dataset_plot_or$p_bonf_gprd >= 0.01 ~ 0.15))
segments(x0=subs$r_squared_4699605, x1=subs$r_squared_4699605+x_label_offset, y0=subs$conservative_or, y1=subs$conservative_or/y_label_offset, col=eqtl_color)
text(subs$r_squared_4699605, subs$conservative_or/y_label_offset, col=eqtl_color, pos=4, labels='rs17327121 (cerebellum eQTL lead SNP)', cex=0.8) 

mtext("D", side=3, cex=1.5, adj = 0, line = 0.5)

write_tsv(all_dataset, 'output/targeted_summary_statistics.tsv.gz')


silence_is_golden = dev.off()
}

## Figure 4 - zoom in on variants around PRNP #### 

if ('figure-4'=='figure-4') {

resx=300
png('display_items/figure-4.png',width=6.5*resx,height=6.5*resx,res=resx)

layout_matrix = matrix(c(1,5,
                         2,6,
                         3,7,
                         4,8), nrow=4, byrow=T)
layout(layout_matrix, heights=c(1,1,1,0.15))


par(mar=c(1,4,3,1))
xlims = c(4570000, 4710000)
ylims = c(pseud0,pseudinf_gprd+10)
xbigs = seq(min(xlims), max(xlims), by=10000)
xats = seq(min(xlims), max(xlims), by=1000)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F,log="y")
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
points(all_dataset_plot_or$postion, all_dataset_plot_or$conservative_or, pch=20, cex=0.5)
abline(h=1, lwd=0.5)
points(all_dataset_plot_or$postion, all_dataset_plot_or$conservative_or, pch=20, cex=0.5,
       col = case_when(all_dataset_plot_or$r_squared_4699605 >= 0.25 ~ linked_color,
                       all_dataset_plot_or$r_squared_4699605 <= -0.25 ~ antilinked_color))
axis(side=2, at=c(pseud0/2,pseudinf_gprd+10), lwd.ticks=0, labels = NA, las=2)
ybigs = c(0.05,0.1,1,10, 20)
ybiglabs = c('0','0.1' ,'1', '10', '20')
axis(side=2, at=ybigs, labels=ybiglabs, las=2)
axis.break(axis=2, breakpos=0.06, style='slash')
yats = rep(1:9, 4) * rep(10^(-2:1), each=9)
yats = yats[7:(length(yats)-7)]
axis(side=2, at=yats, labels=NA, tck=-0.015)
mtext(side=2, line=2.5, text='OR VPSPr', cex=0.8)
subs = all_dataset_plot_or %>% filter(postion == m129v_pos)
points(subs$postion, subs$conservative_or, pch=1, col=m129v_color, lwd=2)
text(subs$postion, subs$conservative_or, pos=3, col=m129v_color, labels='129V')
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xats, tck=-0.02, labels=NA)

mtext("A", side=3, cex=1.5, adj = -0.2, line = 0.5)



signed = function(x) ifelse(x > 0, paste0('+',x), x)
par(mar=c(1,4,1,1))
ylims = c(-0.3, 0.3)
ybigs = seq(min(ylims), max(ylims), by=0.1)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=ybigs, labels=signed(round(ybigs,1)), las=2)
mtext(side=2, line=3.0, text='tibial nerve slope', cex=0.8)
points(all_dataset_plot_or$postion, all_dataset_plot_or$slope_tibial_nerve, pch=20, cex=0.25, col='black')
abline(h=0)
subs = all_dataset_plot_or %>% filter(postion == tibial_nerve_eqtl_lead_pos)
points(subs$postion, subs$slope_tibial_nerve, pch=1, col=eqtl_color, lwd=2)
#text(subs$postion, subs$slope_tibial_nerve, pos=4, col=eqtl_color, labels='rs6052749 (tibial nerve eQTL lead SNP)')


par(xpd=T)

text(subs$postion, subs$slope_tibial_nerve, 
     pos = 4, col = eqtl_color, 
     labels = "rs6052749 (tibial", cex = 0.8)
text(subs$postion, subs$slope_tibial_nerve - 0.03, 
     pos = 4, col = eqtl_color, 
     labels = "nerve eQTL lead SNP)", cex = 0.8)
par(xpd=F)



subs = all_dataset_plot_or %>% filter(postion == m129v_pos)
points(subs$postion, subs$slope_tibial_nerve, pch=1, col=m129v_color, lwd=2)
text(subs$postion, subs$slope_tibial_nerve, pos=3, col=m129v_color, labels='129V')
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xats, tck=-0.02, labels=NA)
mtext("C", side=3, cex=1.5, adj = -0.2, line = 0.5)



par(mar=c(1.5,4,1,1))
ylims = c(-0.8, 0.8)
ybigs = seq(min(ylims), max(ylims), by=0.1)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=ybigs, labels=signed(ybigs), las=2)
mtext(side=2, line=3.0, text='cerebellum slope', cex=0.8)
points(all_dataset_plot_or$postion, all_dataset_plot_or$slope_cerebellum, pch=20, cex=0.25, col='black')
abline(h=0)
subs = all_dataset_plot_or %>% filter(postion == cerebellum_eqtl_lead_pos)
points(subs$postion, subs$slope_cerebellum, pch=1, col=eqtl_color, lwd=2)
#text(subs$postion, subs$slope_cerebellum, pos=4, col=eqtl_color, labels='rs17327121 (cerebellum eQTL lead SNP)')



text(subs$postion, subs$slope_cerebellum, 
     pos = 4, col = eqtl_color, 
     labels = "rs17327121", cex = 0.8)
text(subs$postion, subs$slope_cerebellum - 0.08, 
     pos = 4, col = eqtl_color, 
     labels = "(cerebellum eQTL lead SNP)", cex = 0.8)




subs = all_dataset_plot_or %>% filter(postion == m129v_pos)
points(subs$postion, subs$slope_cerebellum, pch=1, col=m129v_color, lwd=2)
text(subs$postion, subs$slope_cerebellum, pos=3, col=m129v_color, labels='129V')
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, line=-0.5, lwd=0, labels=paste0(formatC(xbigs/1e6,format='f',digits=2), 'M'))
axis(side=1, at=xats, tck=-0.02, labels=NA)
mtext("E", side=3, cex=1.5, adj = -0.2, line = 0.5)



par(mar=c(0.5,4,0,1))

ylims = c(-2, 6)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
segments(x0=landmarks$exon1_start, x1=landmarks$exon3_end, y0=landmarks$y, lwd=intron_lwd)
rect(xleft=landmarks$exon1_start, xright=landmarks$exon1_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=default_fill, border=NA)
rect(xleft=landmarks$exon3_start, xright=landmarks$exon3_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=default_fill, border=NA)
rect(xleft=landmarks$cds_start, xright=landmarks$cds_end, ybottom=landmarks$y-cds_height/2, ytop=landmarks$y+cds_height/2, col=default_fill, border=NA)
# gene symbols
text(x=c(4694022,4725184), y=landmarks$y-3, labels=toupper(landmarks$gene), cex=.7, font=3)
text(x=min(xlims), y=landmarks$y[1], pos=4, labels='chr20 position', cex=0.9)




par(mar=c(3,4,3,1))
xlims = c(0.7, 1.3)
ylims = c(0.05, 20)
xbigs = c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F, log='xy')
ybigs = c(0.05,0.1,1,10, 20)
ybiglabs = c('0','0.1' ,'1', '10', '20')
axis(side=2, at=ybigs, labels=ybiglabs, las=2)
axis.break(axis=2, breakpos=0.06, style='slash')
yats = rep(1:9, 4) * rep(10^(-2:1), each=9)
yats = yats[7:(length(yats)-7)]
axis(side=2, at=yats, labels=NA, tck=-0.015)
mtext(side=1, line=2.5, text='OR sCJD GWAS', cex=0.8)
mtext(side=2, line=2.5, text='OR VPSPr', cex=0.8)
axis(side=1, at=xbigs)
abline(h=1)
abline(v=1)
points(all_dataset_plot_or$or_scjd_gwas, all_dataset_plot_or$conservative_or,  pch=20, cex=0.5,
       col = case_when(all_dataset_plot_or$r_squared_4699605 >= 0.25 ~ linked_color,
                       all_dataset_plot_or$r_squared_4699605 <= -0.25 ~ antilinked_color,
                       TRUE ~ default_color))
subs = all_dataset_plot_or %>% filter(postion == m129v_pos)
par(xpd=T)
points(subs$or_scjd_gwas, subs$conservative_or, pch=1, col=m129v_color, lwd=2)
text(subs$or_scjd_gwas, subs$conservative_or, pos=3, col=m129v_color, labels='129V')
par(xpd=F)
mtext("B", side=3, cex=1.5, adj = 0, line = 0.5)


par(mar=c(3,4,3,1))
subs = all_dataset_plot_or %>%
  #filter(r_squared_4699605 > -0.25 & r_squared_4699605 < 0.25) %>%
  filter(postion < 4640000)
ylims = c(0.05, 20)
xlims = c(-0.3, 0.3)
xbigs = seq(min(xlims), max(xlims), by=0.1)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F, log='y')
points(subs$slope_tibial_nerve, subs$conservative_or, pch=20, cex=0.5)
abline(v=0)
abline(h=1)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=1, at=xbigs, labels=signed(round(xbigs,1)))
mtext(side=1, line=2.5, text='tibial nerve slope', cex=0.8)
axis(side=2, at=c(pseud0/2,pseudinf_gprd+10), lwd.ticks=0, labels = NA, las=2)
ybigs = c(0.05,0.1,1,10, 20)
ybiglabs = c('0','0.1' ,'1', '10', '20')
axis(side=2, at=ybigs, labels=ybiglabs, las=2)
axis.break(axis=2, breakpos=0.06, style='slash')
yats = rep(1:9, 4) * rep(10^(-2:1), each=9)
yats = yats[7:(length(yats)-7)]
axis(side=2, at=yats, labels=NA, tck=-0.015)
mtext(side=2, line=2.5, text='OR VPSPr', cex=0.8)
cor.test(subs$slope_tibial_nerve, subs$conservative_or)
subs = all_dataset_plot_or %>% filter(postion == m129v_pos)
par(xpd=T)
points(subs$slope_tibial_nerve, subs$conservative_or, pch=1, col=m129v_color, lwd=2)
text(subs$slope_tibial_nerve, subs$conservative_or, pos=3, col=m129v_color, labels='129V')
par(xpd=F)
mtext("D", side=3, cex=1.5, adj = 0, line = 0.5)



par(mar=c(3,4,3,1))
subs = all_dataset_plot_or %>%
  #filter(r_squared_4699605 > -0.25 & r_squared_4699605 < 0.25) %>%
  filter(postion < 4640000)
ylims = c(0.05, 20)
xlims = c(-0.8, 0.8)
xbigs = seq(min(xlims), max(xlims), by=0.1)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F, log='y')
points(subs$slope_cerebellum, subs$conservative_or, pch=20, cex=0.5)
abline(v=0)
abline(h=1)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=1, at=xbigs, labels=signed(round(xbigs,1)))
mtext(side=1, line=2.5, text='cerebellum slope', cex=0.8)
axis(side=2, at=c(pseud0/2,pseudinf_gprd+10), lwd.ticks=0, labels = NA, las=2)
ybigs = c(0.05,0.1,1,10, 20)
ybiglabs = c('0','0.1' ,'1', '10', '20')
axis(side=2, at=ybigs, labels=ybiglabs, las=2)
axis.break(axis=2, breakpos=0.06, style='slash')
yats = rep(1:9, 4) * rep(10^(-2:1), each=9)
yats = yats[7:(length(yats)-7)]
axis(side=2, at=yats, labels=NA, tck=-0.015)
mtext(side=2, line=2.5, text='OR VPSPr', cex=0.8)
cor.test(subs$slope_cerebellum, subs$conservative_or)
subs = all_dataset_plot_or %>% filter(postion == m129v_pos)
par(xpd=T)
points(subs$slope_cerebellum, subs$conservative_or, pch=1, col=m129v_color, lwd=2)
text(  subs$slope_cerebellum, subs$conservative_or, pos=3, col=m129v_color, labels='129V')
par(xpd=F)
mtext("F", side=3, cex=1.5, adj = 0, line = 0.5)


silence_is_golden = dev.off()
}


## Exploratory analyses #### 

# 
# 
# 
# 
# points(all_dataset_plot_or$postion, all_dataset_plot_or$conservative_or, pch=20, cex=0.5,
#        col = case_when(all_dataset_plot_or$r_squared_4699605 >= 0.25 ~ linked_color,
#                        all_dataset_plot_or$r_squared_4699605 <= -0.25 ~ antilinked_color,
#                        all_dataset_plot_or$slope_tibial_nerve > 0.1 ~ '#00FFFF',
#                        all_dataset_plot_or$slope_tibial_nerve < -0.1 ~ '#0088FF',
#                        all_dataset_plot_or$slope_cerebellum > 0.1 ~  '#FF00FF',
#                        all_dataset_plot_or$slope_cerebellum < -0.1 ~ '#880088'))
# points(all_dataset_plot_or$postion, all_dataset_plot_or$conservative_or, pch=20, cex=0.5,
#        col = case_when(all_dataset_plot_or$r_squared_4699605 >= 0.25 ~ linked_color,
#                        all_dataset_plot_or$r_squared_4699605 <= -0.25 ~ antilinked_color))


# exploring


all_dataset_plot_or %>%
 filter((conservative_or < 0.2 | conservative_or > 5) & abs(r_squared_4699605) < 0.3) -> lookup_variants
# this yields 2 variants with OR < 0.2:

# https://gnomad.broadinstitute.org/variant/20-4686278-C-G?dataset=gnomad_r4
# "Warning This variant is covered in fewer than 50% of individuals in gnomAD v4.1.0 exomes. Allele frequency estimates may not be reliable."
# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A4686049-4686507
# highly GC-rich region immediately upstream of PRNP promoter

# https://gnomad.broadinstitute.org/variant/20-4692886-G-A?dataset=gnomad_r4
# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A4674297-4711475
# also in a fairly GC-rich region in the middle of the intron, right inbetween two LINE repeats

# and two with OR > 5:
# https://gnomad.broadinstitute.org/variant/20-4630961-C-T?dataset=gnomad_r4
# nothing obviously wrong with it.
# https://gnomad.broadinstitute.org/variant/20-4701686-C-T?dataset=gnomad_r4
# nothing obviously wrong with it.


