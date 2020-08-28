###
##
##            Longevity Returns to Political Office
##
###

sink("analysis-output.txt") 

# Load packages
library("dplyr")
library("ggplot2")
library("purrr")
library("texreg") 
library("data.table") 
library("rdrobust") 
library("stargazer")
library("scales")
library("digest")

set.seed(100)

# Set ggplot2 theme
theme_set(
  theme_grey(base_size = 11.5) %+replace% 
    theme(
      plot.margin = unit(rep(0.5, 4), "cm"), plot.background = element_blank(), panel.background = element_blank(),
      panel.border = element_blank(), legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA), legend.title = element_blank(),
      strip.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_line(linetype = "dotted", colour = "#757575", size = 0.3), panel.grid.minor = element_blank(),
      axis.ticks = element_blank(), axis.line = element_line(color = "#FFFFFF", size = 0.3),
      plot.title = element_text(size = 12, hjust = 0, margin = margin(b = 15)),
      plot.subtitle = element_text(size = 12, hjust = 0, margin = margin(b = 5)),
      plot.caption = element_text(size = 10, colour = "#212121", margin = margin(t = 15)),
      axis.title = element_text(size = 11, face = "plain"), axis.text = element_text(size = 10, face = "plain"),
      legend.text = element_text(size = 10), strip.text = element_text(size = 12, face = "plain")
    )
)

# Functions
## 'stardd': Export 'rdrobust' object to latex
stardd = function(rdd, coef_name = "LATE"){
  # Input: 'rdrobust' object
  # rdd: list object
  
  if( rdd$bws["h", "left"] != rdd$bws["h", "right"]){
    stop("Not constant bandwith")
  }
  
  else{
    
    n = sum(rdd$N)
    n_eff = rdd$Nh[1] + rdd$Nh[2]
    
    # Extract bandwith
    bw = rdd$bws["h", "left"]
    
    # Extract coef + SE
    coef = rdd$coef[3]
    se = rdd$se[3]
    pval = rdd$pv[3]
    
    tr = createTexreg(
      coef.names = coef_name,
      coef = coef,
      se = se, 
      pvalues = pval,
      gof = c(n, n_eff, bw),
      gof.names = c("Observations", "Effective observations", "Bandwidth"),
      gof.decimal = c(FALSE, FALSE, TRUE)
    )
    return(tr)
  }
}

## Export to latex
rdd_to_latex = function(x,x_names = "Election win", ...){
  # x: list of texreg objects
  texreg(x, 
         ...,
         stars = c(0.01, 0.05, 0.1),
         custom.coef.names = x_names,
         custom.note = "Calonico et al. (2014) optimal bandwith and triangular kernel weights in all columns. All models use local linear regression and include the bias correction and robust standard errors of Calonico et al. (2014). %stars.")
}

tidy_rdd = function(rdd, coef_name = "LATE"){
  # Input: 'rdrobust' object
  # rdd: list object
  
  if( rdd$bws["h", "left"] != rdd$bws["h", "right"]){
    stop("Not constant bandwith")
  }
  
  else{
    
    n = rdd$N
    n_eff = rdd$Nh[1] + rdd$Nh[2]
    
    # Extract bandwith
    bw = rdd$bws["h", "left"]
    
    # Extract coef + SE
    coef = rdd$coef[3]
    se = rdd$se[3]
    z = rdd$z[3]
    pval = rdd$pv[3]
    
    x = data.frame(
      coef = coef,
      se = se,
      z = z,
      pval = pval,
      stringsAsFactors = FALSE
    )
    
    names(x) = c("Estimate",
                 "Std. Error",
                 "$Z$ value",
                 "$P$ value")
    
    x
    
  }
}

df_rdd <- fread("longevity.csv")

if (sha1(df_rdd) != "300cc29bbecd2b630016c9bd2c8ef958dcc1b45d"){
  error("Wrong data file loaded or data has been changed!") }

df_m <- df_rdd %>% 
  filter(year >= 1945, living_day_imp_post > 0) 

# Numerical statements
## Still alive as of September 2019"
cat("Numbers reported in manuscript:\n")
cat(paste0("Still living: ", NROW(df_m[df_rdd$year >= 1945 & df_rdd$living == "yes" & !is.na(df_rdd$living),]), "\n\n"))

## Candidate-year observations
cat(paste0("Candidate-year observations: ", NROW(df_m), "\n\n"))

# Number of elections:
n_elections <- df_m %>% count(area, year) %>% NROW()
cat(paste0("Number of unique elections: ", n_elections, "\n\n"))

# Number of D, R, Third
cat(paste0("Number of Democrats: ", table(df_m$party)[1], "\n"))
cat(paste0("Number of Republicans: ", table(df_m$party)[2], "\n"))
cat(paste0("Number of Third party: ", table(df_m$party)[3], "\n\n"))

# Number of times running for office
cat("Number of times running for office:\n")
df_rdd %>%
  filter(year >= 1945) %>%
  count(cand_last, cand_first) %>%
  count(n) %>%
  mutate(prop = nn / sum(nn))

### Table 1
main_1 <- rdrobust(y = df_m$living_day_imp_post, x = df_m$margin_pct_1)

main_2 <- rdrobust(y = df_m$living_day_imp_post, x = df_m$margin_pct_1,
                  covs = cbind(df_m$reg_northeast, df_m$reg_south, df_m$reg_west))

main_3 <- rdrobust(y = df_m$living_day_imp_post, x = df_m$margin_pct_1,
                  covs = cbind(df_m$female, df_m$ex, df_m$republican, df_m$third))

main_4 <- rdrobust(y = df_m$living_day_imp_post, x = df_m$margin_pct_1,
                  covs = cbind(df_m$female, df_m$ex, df_m$republican, df_m$third, 
                               df_m$reg_northeast, df_m$reg_south, df_m$reg_west)) 


cat("\nTable 1:\n")
list(main_1, main_2, main_3, main_4) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "table1.tex"

h <- main_1$bws["h", "left"]

# Plot main result
x <- rdplot(y = df_m$living_day_imp_post,
           x = df_m$margin_pct_1,
           h = h,
           subset = df_m$margin_pct_1<=h & 
             df_m$margin_pct_1>=-h,
           binselect = "esmv", kernel="triangular", p=1,
           title = "Main result",
           y.label = "Days alive after election",
           x.label = "Vote Share")

df_plot <- x$vars_bins

df_plot_2 <- data.frame(numbers = seq(min(df_plot$rdplot_mean_x), 
                max(df_plot$rdplot_mean_x), 
                by = 0.01)
)

df_plot_2$line <- c(x$coef[1, 1] + x$coef[2, 1]*df_plot_2$numbers[df_plot_2$numbers < 0],
                    x$coef[1, 2] + x$coef[2, 2]*df_plot_2$numbers[df_plot_2$numbers > 0])

p <- ggplot() +
  geom_line(color = "black", size = 1, data = filter(df_plot_2, numbers <= 0), aes(x = numbers, y = line)) +
  geom_line(color = "black", size = 1, data = filter(df_plot_2, numbers >= 0), aes(x = numbers, y = line)) +
  coord_cartesian(ylim = c(0, 20000)) +
  scale_y_continuous(breaks = seq(0, 15000, by = 5000), labels = scales::comma, expand = c(0, 0)) +
  scale_x_continuous(breaks = c(seq(-10, 10, by = 2.5)), labels = c("-10", "-7.5", "-5", "-2.5", "0\n(Cutoff)", "2.5", "5", "7.5", "10")) +
  geom_vline(xintercept = 0, color = "grey40",size = 1/2, linetype = "dashed") +
  labs(x = "Margin (pct)", y = "Days alive after election") +
  geom_point(alpha = .85, size = 2.5, shape = 21, fill = "grey98", color = "black", data = df_plot, aes(x = rdplot_mean_x, y = rdplot_mean_y))

ggsave(plot = p, file = "figure1.pdf", height = 4, width = 8)

# Show figure 
p

## Appendix

## Appendix Table 1
df_desc = df_m %>% 
  select(
    living_day_imp_post, living_day_imp_pre, living_day_pre, ex, female, democrat, 
    republican, pc_inc_ann, pop_annual, tot_expenditure, reg_south, reg_west, 
    reg_midwest, reg_northeast) %>% 
  mutate(ex = ex*365) %>% 
  as.data.frame()

names(df_desc) = c(
  "Days alive after election",
  "Days alive before election (imputed)",
  "Days alive before election (not imputed)",
  "Life expectancy",
  "Female", "Democrat", "Republican",
  "Per capita income",
  "Population",
  "Total expenditure",
  "Census region: South",
  "Census region: West",
  "Census region: Northeast",
  "Census region: Midwest")

cat("\nAppendix Table 1:\n")
stargazer(df_desc, digits = 2, 
          summary = TRUE, 
          rownames = TRUE,
          # To save file, add: out = "appendix-table1.tex", 
          font.size = "scriptsize",
          label = "a_summary",
          title = "Summary statistics")

## Appendix Table 2
a_1 = rdrobust(y = df_m$pc_inc_ann, x = df_m$margin_pct_1)
a_2 = rdrobust(y = df_m$pop_annual, x = df_m$margin_pct_1)
a_3 = rdrobust(y = df_m$tot_expenditure, x = df_m$margin_pct_1)
a_4 = rdrobust(y = df_m$democrat, x = df_m$margin_pct_1)
a_5 = rdrobust(y = df_m$republican, x = df_m$margin_pct_1)
a_6 = rdrobust(y = df_m$living_day_imp_pre, x = df_m$margin_pct_1)
a_7 = rdrobust(y = df_m$living_day_pre, x = df_m$margin_pct_1)
a_8 = rdrobust(y = df_m$female, x = df_m$margin_pct_1)
a_9 = rdrobust(y = df_m$reg_south, x = df_m$margin_pct_1)
a_10 = rdrobust(y = df_m$reg_west, x = df_m$margin_pct_1)
a_11 = rdrobust(y = df_m$reg_northeast, x = df_m$margin_pct_1)
a_12 = rdrobust(y = df_m$reg_midwest, x = df_m$margin_pct_1)
a_13 = rdrobust(y = df_m$ex*365, x = df_m$margin_pct_1)

x = list(a_4, a_5, a_8, a_6, a_7, a_13, a_1, a_2, a_3, a_9, a_10, a_11, a_12) %>% 
  map_df(tidy_rdd) 

rownames(x) = c("Democrat", "Republican", "Female", "Days alive before election (imputed)",
                "Days alive before election (not imputed)", "Life expectancy", "Per capita income",
                "Population", "Total expenditure", "Census region: South (dummy)",
                "Census region: West (dummy)", "Census region: Northeast (dummy)",
                "Census region: Midwest (dummy)")

cat("\nAppendix Table 2:\n")
stargazer(x, summary = FALSE, rownames = TRUE) # To save file, add: out = "appendix-table2.tex"


## Appendix Table 3
main_1_noim = rdrobust(y = df_m$living_day_post, x = df_m$margin_pct_1)

main_2_noim = rdrobust(y = df_m$living_day_post, x = df_m$margin_pct_1,
                  covs = cbind(df_m$reg_northeast, df_m$reg_south, df_m$reg_west))

main_3_noim = rdrobust(y = df_m$living_day_post, x = df_m$margin_pct_1,
                  covs = cbind(df_m$female, df_m$ex, df_m$republican, df_m$third))

main_4_noim = rdrobust(y = df_m$living_day_post, x = df_m$margin_pct_1,
                  covs = cbind(
                    df_m$female,
                    df_m$ex,
                    df_m$republican,
                    df_m$third,
                    df_m$reg_northeast,
                    df_m$reg_south,
                    df_m$reg_west)) 

cat("\nAppendix Table 3:\n")
list(main_1_noim, main_2_noim, main_3_noim, main_4_noim) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "appendix-table3.tex"

## Appendix Figure 1
coef_df_1 <- 1:20 %>% 
  map(function(x){
    rdrobust(
      y=df_m$living_day_imp_post,
      x=df_m$margin_pct_1,
      h = x, rho = .5)
  }) %>% 
  map_df(tidy_rdd) %>% 
  mutate(bw = 1:20) %>% 
  mutate(
    ind = "Model 1"
  )

coef_df_2 <- 1:20 %>% 
  map(function(x){
    rdrobust(
      y=df_m$living_day_imp_post,
      x=df_m$margin_pct_1,
      h = x, rho = .5,
      covs = cbind(
        df_m$reg_northeast,
        df_m$reg_south,
        df_m$reg_west))
  }) %>% 
  map_df(tidy_rdd) %>% 
  mutate(bw = 1:20) %>% 
  mutate(
    ind = "Model 2"
  )

coef_df_3 <- 2:20 %>% 
  map(function(x){
    rdrobust(
      y=df_m$living_day_imp_post,
      x=df_m$margin_pct_1,
      h = x, rho = .5,
      covs = cbind(
        df_m$female,
        df_m$ex,
        df_m$republican,
        df_m$third))
  }) %>% 
  map_df(tidy_rdd) %>% 
  mutate(bw = 2:20) %>% 
  mutate(
    ind = "Model 3"
  )

coef_df_4 <- 2:20 %>% 
  map(function(x){
    rdrobust(
      y=df_m$living_day_imp_post,
      x=df_m$margin_pct_1,
      h = x, rho = .5,
      covs = cbind(
        df_m$female,
        df_m$ex,
        df_m$republican,
        df_m$third,
        df_m$reg_northeast,
        df_m$reg_south,
        df_m$reg_west)) 
  }) %>% 
  map_df(tidy_rdd) %>% 
  mutate(bw = 2:20) %>% 
  mutate(
    ind = "Model 4"
  )

coef_df <- coef_df_1 %>% 
  bind_rows(coef_df_2) %>% 
  bind_rows(coef_df_3) %>% 
  bind_rows(coef_df_4)

p <- ggplot(
  filter(coef_df, ind %in% c("Model 1", "Model 4")), 
  aes(x = bw, y = Estimate, fill = ind))
p <- p + 
  geom_hline(
    yintercept = 0, size = 0.8, color = "black", linetype = "dashed") +
  geom_linerange(
    position=position_dodge(.6),
    aes(ymin = Estimate - 1.96*`Std. Error`, 
        ymax = Estimate + 1.96*`Std. Error`),
    size = .65) +
  geom_point(
    shape = 21, color = "black",
    position=position_dodge(.6),
    size = 2.75) +
  scale_fill_manual(values = c("black", "grey99"), name = NULL) +
  labs(x = "Bandwidth", y = "Estimate") +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, 6000, by = 3000),
    labels = c("0", "3,000", "6,000")) +
  expand_limits(y = c(-2000, 9000)) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 11, 12, 14, 16, 18, 20),
                     labels = c("2", "4", "6", "8", "10", 
                                "11\n(Optimal)\n", "12", "14", "16",
                                "18", "20")) +
  theme(
    legend.position=c(.925, .9),
    legend.text = element_text(size = 9)
  )

ggsave(plot = p, file = "appendix-figure1.pdf", height = 4, width = 8)


## Appendix Table 4
df_out = df_m %>% 
  filter(living_day_imp_post > quantile(living_day_imp_post, .02, na.rm = TRUE), 
         living_day_imp_post < quantile(living_day_imp_post, .98, na.rm = TRUE))

out_1 = rdrobust(y = df_out$living_day_imp_post, x = df_out$margin_pct_1)

out_2 = rdrobust(y = df_out$living_day_imp_post, 
                 x = df_out$margin_pct_1,
                 covs = cbind(df_out$reg_northeast, df_out$reg_south, df_out$reg_west))

out_3 = rdrobust(y = df_out$living_day_imp_post, 
                 x = df_out$margin_pct_1,
                 covs = cbind(df_out$female, df_out$ex, df_out$republican, df_out$third))

out_4 = rdrobust(y = df_out$living_day_imp_post, 
                 x = df_out$margin_pct_1,
                 covs = cbind(df_out$female, df_out$ex, df_out$republican, df_out$third,
                   df_out$reg_northeast, df_out$reg_south, df_out$reg_west)) 


cat("\nAppendix Table 4:\n")
list(out_1, out_2, out_3, out_4) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "appendix-table4.tex"

## Appendix Table 5
df_m_chall <- df_m %>% 
  filter(status == "challenger")

chall_1 = rdrobust(y = df_m_chall$living_day_imp_post, 
                   x = df_m_chall$margin_pct_1)

chall_2 = rdrobust(y = df_m_chall$living_day_imp_post, 
                   x = df_m_chall$margin_pct_1,
                   covs = cbind(df_m_chall$reg_northeast, df_m_chall$reg_schallh, df_m_chall$reg_west))

chall_3 = rdrobust(y = df_m_chall$living_day_imp_post, 
                   x = df_m_chall$margin_pct_1,
                   covs = cbind(df_m_chall$female, df_m_chall$living_day_imp_pre, df_m_chall$republican, df_m_chall$third))

chall_4 = rdrobust(y = df_m_chall$living_day_imp_post, 
                   x = df_m_chall$margin_pct_1,
                   covs = cbind(df_m_chall$female, df_m_chall$living_day_imp_pre, df_m_chall$republican,
                     df_m_chall$third, df_m_chall$reg_northeast, df_m_chall$reg_schallh, df_m_chall$reg_west)) 

df_m_chall_out = df_m %>% 
  filter(living_day_imp_post > quantile(living_day_imp_post, .02, na.rm = TRUE),
         living_day_imp_post < quantile(living_day_imp_post, .98, na.rm = TRUE),
         status == "challenger") 

chall_out_1 = rdrobust(y = df_m_chall_out$living_day_imp_post, x = df_m_chall_out$margin_pct_1)

chall_out_2 = rdrobust(y = df_m_chall_out$living_day_imp_post, 
                       x = df_m_chall_out$margin_pct_1,
                       covs = cbind(df_m_chall_out$reg_northeast, df_m_chall_out$reg_south, df_m_chall_out$reg_west))

chall_out_3 = rdrobust(y = df_m_chall_out$living_day_imp_post, 
                       x = df_m_chall_out$margin_pct_1,
                       covs = cbind(
                         df_m_chall_out$female,
                         df_m_chall_out$living_day_imp_pre,
                         df_m_chall_out$republican,
                         df_m_chall_out$third))

chall_out_4 = rdrobust(y = df_m_chall_out$living_day_imp_post, 
                       x = df_m_chall_out$margin_pct_1,
                       covs = cbind(
                         df_m_chall_out$female,
                         df_m_chall_out$living_day_imp_pre,
                         df_m_chall_out$republican,
                         df_m_chall_out$third,
                         df_m_chall_out$reg_northeast,
                         df_m_chall_out$reg_south,
                         df_m_chall_out$reg_west)) 

cat("\nAppendix Table 5:\n")
list(chall_1, chall_2, chall_3, chall_4) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "appendix-table5a.tex"

list(chall_out_1, chall_out_2, chall_out_3, chall_out_4) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "appendix-table5b.tex"

## Appendix Table 6
main_1_2nd = rdrobust(y = df_m$living_day_imp_post, x = df_m$margin_pct_1, p = 2)

main_2_2nd = rdrobust(y = df_m$living_day_imp_post, x = df_m$margin_pct_1, p = 2,
                      covs = cbind(
                        df_m$reg_northeast,
                        df_m$reg_south,
                        df_m$reg_west))

main_3_2nd = rdrobust(y = df_m$living_day_imp_post, x = df_m$margin_pct_1, p = 2,
                      covs = cbind(
                        df_m$female,
                        df_m$ex,
                        df_m$republican,
                        df_m$third))

main_4_2nd = rdrobust(y = df_m$living_day_imp_post, x = df_m$margin_pct_1, p = 2,
                      covs = cbind(
                        df_m$female,
                        df_m$ex,
                        df_m$republican,
                        df_m$third,
                        df_m$reg_northeast,
                        df_m$reg_south,
                        df_m$reg_west)) 


cat("\nAppendix Table 6:\n")
list(main_1_2nd, main_2_2nd, main_3_2nd, main_4_2nd) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "appendix-table6.tex"

## Appendix Figure 2
df_m_neg = df_m %>% 
  filter(margin_pct_1 < 0)

out_1 = rdrobust(y = df_m_neg$living_day_imp_post, 
                 x = df_m_neg$margin_pct_1,
                 c = -1)

out_2 = rdrobust(y = df_m_neg$living_day_imp_post, 
                 x = df_m_neg$margin_pct_1,
                 c = -5)

out_3 = rdrobust(y = df_m_neg$living_day_imp_post, 
                 x = df_m_neg$margin_pct_1,
                 c = -10)

df_m_pos = df_m %>% 
  filter(margin_pct_1 > 0)

out_1 = rdrobust(y = df_m_pos$living_day_imp_post, 
                 x = df_m_pos$margin_pct_1,
                 c = 1)

out_2 = rdrobust(y = df_m_pos$living_day_imp_post, 
                 x = df_m_pos$margin_pct_1,
                 c = 5)

out_3 = rdrobust(y = df_m_pos$living_day_imp_post, 
                 x = df_m_pos$margin_pct_1,
                 c = 10)


rdd_function = function(c){
  rdrobust(y = df_m$living_day_imp_post, 
           x = df_m$margin_pct_1,
           c = c, h = 10)
}

cutoff = seq(-25, 25, by = 1) 
df_plot = cutoff %>% 
  map(rdd_function) %>% 
  map_df(tidy_rdd)  %>% 
  mutate(c = cutoff)

p = ggplot(
  df_plot, 
  aes(x = c, y = Estimate))
p = p +
  geom_hline(yintercept = 0, size = 0.8, linetype="dashed") +
  geom_ribbon(aes(ymin = Estimate - 1.96 * `Std. Error`,
                  ymax = Estimate + 1.96 * `Std. Error`),
              alpha = .45) +
  geom_line() +
  labs(y = "Estimate", x = "Cutoff") +
  scale_y_continuous(labels = scales::comma) 

ggsave(plot = p, file = "appendix-figure2.pdf", height = 4, width = 8)

## Appendix Table 7
df_old = df_rdd %>% 
  filter(living == "no") %>% 
  filter(year >= 1945 & year < 1970) %>% 
  filter(living_day_imp_post > 0) 

main_old_1 = rdrobust(y = df_old$living_day_imp_post, 
                      x = df_old$margin_pct_1)

main_old_2 = rdrobust(y = df_old$living_day_imp_post, 
                      x = df_old$margin_pct_1,
                      covs = cbind(
                        df_old$reg_northeast,
                        df_old$reg_south,
                        df_old$reg_west))

main_old_3 = rdrobust(y = df_old$living_day_imp_post, 
                      x = df_old$margin_pct_1,
                      covs = cbind(
                        df_old$ex,
                        df_old$republican
                      ))


main_old_4 = rdrobust(y = df_old$living_day_imp_post, 
                      x = df_old$margin_pct_1,
                      covs = cbind(
                        df_old$ex,
                        df_old$republican,
                        df_old$reg_northeast,
                        df_old$reg_south,
                        df_old$reg_west)) 


cat("\nAppendix Table 7:\n")
list(main_old_1, main_old_2, main_old_3, main_old_4) %>% 
  map(stardd) %>% 
  rdd_to_latex() # to save file, add: file = "appendix-table7.tex"

## Appendix Table 8
df_missing = df_rdd %>% 
  filter(year >= 1945) %>% 
  mutate(death_missing = ifelse(is.na(death_date_imp), 1, 0)) 

summary(rdrobust(y = df_missing$death_missing, 
                 x = df_missing$margin_pct_1))

df_nomissing <- df_rdd %>%
  filter(year >= 1945) %>%
  group_by(area) %>%
  filter(!is.na(death_date_imp)) %>%
  count(year) %>%
  mutate(missdeath = ifelse(n == 2, 0, 1)) %>%
  select(-n)

df_nomissdata <- left_join(df_rdd, df_nomissing, by = c("area", "year")) %>%
  filter(missdeath == 0)

main_nomiss_1 = rdrobust(y = df_nomissdata$living_day_imp_post, 
                         x = df_nomissdata$margin_pct_1)

main_nomiss_2 = rdrobust(y = df_nomissdata$living_day_imp_post, 
                         x = df_nomissdata$margin_pct_1,
                         covs = cbind(
                           df_nomissdata$reg_northeast,
                           df_nomissdata$reg_south,
                           df_nomissdata$reg_west))

main_nomiss_3 = rdrobust(y = df_nomissdata$living_day_imp_post, 
                         x = df_nomissdata$margin_pct_1,
                         covs = cbind(
                           df_nomissdata$female,
                           df_nomissdata$ex,
                           df_nomissdata$republican))


main_nomiss_4 = rdrobust(y = df_nomissdata$living_day_imp_post, 
                         x = df_nomissdata$margin_pct_1,
                         covs = cbind(
                           df_nomissdata$female,
                           df_nomissdata$ex,
                           df_nomissdata$republican,
                           df_nomissdata$reg_northeast,
                           df_nomissdata$reg_south,
                           df_nomissdata$reg_west)) 

cat("\nAppendix Table 8:\n")
list(main_nomiss_1, main_nomiss_2, main_nomiss_3, main_nomiss_4) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "appendix-table8.tex"

## Appendix Table 9
df_m_r = df_m %>% 
  filter(party == "r")

df_m_d = df_m %>% 
  filter(party == "d")

out_1_r = rdrobust(y = df_m_r$living_day_imp_post, 
                   x = df_m_r$margin_pct_1)

out_2_r = rdrobust(y = df_m_r$living_day_imp_post, 
                   x = df_m_r$margin_pct_1,
                   covs = cbind(
                     df_m_r$reg_northeast,
                     df_m_r$reg_south,
                     df_m_r$reg_west))

out_3_r = rdrobust(y = df_m_r$living_day_imp_post, 
                   x = df_m_r$margin_pct_1,
                   covs = cbind(
                     df_m_r$living_day_imp_pre))

out_4_r = rdrobust(y = df_m_r$living_day_imp_post, 
                   x = df_m_r$margin_pct_1,
                   covs = cbind(
                     df_m_r$living_day_imp_pre,
                     df_m_r$reg_northeast,
                     df_m_r$reg_south,
                     df_m_r$reg_west)) 


out_1_d = rdrobust(y = df_m_d$living_day_imp_post, 
                   x = df_m_d$margin_pct_1)

out_2_d = rdrobust(y = df_m_d$living_day_imp_post, 
                   x = df_m_d$margin_pct_1,
                   covs = cbind(
                     df_m_d$reg_northeast,
                     df_m_d$reg_south,
                     df_m_d$reg_west))

out_3_d = rdrobust(y = df_m_d$living_day_imp_post, 
                   x = df_m_d$margin_pct_1,
                   covs = cbind(
                     df_m_d$living_day_imp_pre))

out_4_d = rdrobust(y = df_m_d$living_day_imp_post, 
                   x = df_m_d$margin_pct_1,
                   covs = cbind(
                     df_m_d$living_day_imp_pre,
                     df_m_d$reg_northeast,
                     df_m_d$reg_south,
                     df_m_d$reg_west)) 

cat("\nAppendix Table 9:\n")
list(out_1_r, out_2_r, out_3_r, out_4_r) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "appendix-table9a.tex"

list(out_1_d, out_2_d, out_3_d, out_4_d) %>% 
  map(stardd) %>% 
  rdd_to_latex() # To save file, add: file = "appendix-table9b.tex"

cat("\nSession Info:\n")
sessionInfo()

sink() 

# Save README.txt file
sink("README.txt") 
cat("\nThe results produced with longevity.R and the data in longevity.csv with this session in R:\n\n")
sessionInfo()
sink() 