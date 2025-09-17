# libraries ---------------------------------------------------------------
library(tidyverse)


# read in the data --------------------------------------------------------
# dataset senmayo
df_senmayo <- read_tsv("../../out/table/revision/121_SENMAYO_prop_table.tsv")  

# confirm this has the zeroes
df_senmayo %>%
  filter(prop == 0)

# dataset SIT
df_SIT <- read_tsv("../../out/table/revision/124_SIT_prop_table.tsv")  

# confirm this has the zeroes
df_SIT %>%
  filter(prop == 0)

# wrangling ---------------------------------------------------------------
# are there missing samples
full_join(df_senmayo %>%
            group_by(orig.ident,origin,pathology_class.x) %>%
            summarise(n = n()),
          df_SIT %>%
            group_by(orig.ident,origin,pathology_class.x) %>%
            summarise(n = n()),by = c("orig.ident","origin","pathology_class.x"),suffix = c(".senmayo",".SIT")) %>%
  mutate(delta = n.senmayo - n.SIT) %>%
  filter(delta != 0)

# confirm that the total number of cell is the same for both 
full_join(read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000_expertAnno/enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise.tsv") %>%
            filter(signature == "senmayo") %>%
            group_by(orig.ident,expertAnno.l1,pathology_class) %>%
            summarise(tot = sum(n)),
          read_tsv("../../out/table/ManualClean/enumeration_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno.tsv") %>%
            group_by(orig.ident,expertAnno.l1,pathology_class) %>%
            summarise(tot = sum(n)),by = c("orig.ident","expertAnno.l1","pathology_class"),suffix = c(".senmayo",".SIT")) %>%
  mutate(delta = tot.senmayo - tot.SIT) %>%
  filter(delta != 0)


#
df_plot <- full_join(df_senmayo %>%
                       select(orig.ident,origin,pathology_class.x,expertAnno.l1,n,tot,prop),
                     df_SIT %>%
                       select(orig.ident,origin,pathology_class.x,expertAnno.l1,n,tot,prop),by = c("orig.ident","origin","pathology_class.x","expertAnno.l1"),suffix = c(".senmayo",".SIT"))

# save the plot
df_plot %>%
  write_tsv("../../out/table/124_correlation_SIT_SENMAYO.tsv")

# plot the correlation
# notice that in this case I am discarding the samples for which there is a estimate of 0 senescence
df_plot %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  facet_wrap(~expertAnno.l1,scales = "free")+theme(strip.background = element_blank())
ggsave(paste0("../../out/image/124_correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_split.pdf"),width = 9,height =8)

# plot a globa correlation
df_plot %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/image/124_correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_global.pdf"),width = 5,height =5)

# remove the zeros
df_plot %>%
  filter(prop.SIT!=0,
         prop.senmayo!=0) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/image/124_correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_global_noZero.pdf"),width = 5,height =5)
