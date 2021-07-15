#### R code for generating awesome plot
library(ggplot2)
library(ggthemes)
library(reshape2)
#######################################################################
## Dot plot each sample with each lane
df <- data.frame(Treatment= rep(c("A", "B", "C"), each = 100),
                 LocA=sample(1:100, 50), 
                 LocB=sample(1:150, 50), 
                 LocC=sample(1:100, 50),
                 LocD=sample(1:300, 50),
                 LocE=sample(1:500, 50),
                 locD=sample(1:350, 50))
df_melted <- melt(df, id.vars=c("Treatment"))
colnames(df_melted)[2] <- "Samples"
ggplot(data=df_melted, aes(x=Samples, y=value, color=Treatment)) + 
  geom_point(position=position_dodge(width=0.3))+
  labs(title="Dot plot", 
       subtitle="Drug treatment response",
       caption="Drug : mg", ## This will like notes/citations
       x="Drugs",
       y="Doses")

#######################################################################
## Here I am deciding size of dots on basis of value and color based on Treatment
## We removed major and minor grid background white (theme_bw())
df <- data.frame(Treatment= rep(c("A", "B", "C"), each = 50),
                 LocA=sample(1:100, 50), 
                 LocB=sample(1:150, 50), 
                 LocC=sample(1:100, 50),
                 LocD=sample(1:300, 50),
                 LocE=sample(1:500, 50),
                 LocE=sample(1:350, 50))
df_melted <- melt(df, id.vars=c("Treatment"))
colnames(df_melted)[2] <- "Samples"
ggplot(data=df_melted, aes(x=Samples, y=value, color=Treatment)) + 
  geom_point(aes(col=Treatment, size=value),position=position_dodge(width=0.3))+
  labs(title="Dot plot", 
       subtitle="Drug treatment response",
       caption="Drug : mg",
       x="Drugs",
       y="Doses") +theme_bw()+ # either I have to use this line or plot.background = element_blank(). Both are same
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
######################################################################
## Same dot plot, but this one have violin plot and axix is flipped
## Both title and subtitles are center aligned (by default it's left aligned)
ggplot(data=df_melted, aes(x=Samples, y=value, color=Treatment)) + 
  geom_violin(alpha=0.5, color="gray")+
  geom_point(aes(col=Treatment, size=value),position=position_dodge(width=0.3))+coord_flip()+
  labs(title="Dot plot", 
       subtitle="Drug treatment response",
       caption="Drug : mg",
       x="Drugs",
       y="Doses")+theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
####################################################################
## Same plot with violin plot, but jittereed point (means point are close to each other)
## It's look more organized
ggplot(data=df_melted, aes(x=Samples, y=value, color=Treatment)) + 
  geom_violin(alpha=0.5, color="gray")+ 
  geom_jitter(alpha=0.5, aes(col=Treatment, size=value),position=position_jitter(width=0.1))+coord_flip()+
  labs(title="Dot plot", 
       subtitle="Drug treatment response",
       caption="Drug : mg",
       x="Drugs",
       y="Doses")+theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
######################################################################
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
ggplot(data=df_melted, aes(x=Samples, y=value, color=Treatment)) + 
  geom_violin(alpha=0.5, color="gray")+ 
  #geom_boxplot(width=0.2)+
  #stat_summary(fun.data=data_summary)+
  geom_jitter(alpha=0.5, aes(col=Treatment, size=value),position=position_jitter(width=0.1))+coord_flip()+
  labs(title="Dot plot", 
       subtitle="Drug treatment response",
       caption="Drug : mg",
       x="Drugs",
       y="Doses")+theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #geom_boxplot(width=0.2)+ ## This will add boxplot in violin plot
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
######################################################################
## same figure with box plot
ggplot(data=df_melted, aes(x=Samples, y=value, color=Treatment)) + 
  geom_boxplot(alpha=0.5, color="gray")+ 
  geom_jitter(alpha=0.5, aes(col=Treatment, size=value),position=position_jitter(width=0.1))+coord_flip()+
  #geom_smooth(method="loess", se=F) + ## This will be use for correlation line
  labs(title="Dot plot", 
       subtitle="Drug treatment response",
       caption="Drug : mg",
       x="Drugs",
       y="Doses")+theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
######################################################################
#same plot as over it but notch=TRUE in boxplot
ggplot(data=df_melted, aes(x=Samples, y=value, color=Treatment)) + 
  geom_boxplot(alpha=0.5, color="gray", notch = TRUE)+ 
  geom_jitter(alpha=0.5, aes(col=Treatment, size=value),position=position_jitter(width=0.1))+coord_flip()+
  #geom_smooth(method="loess", se=F) + ## This will be use for correlation line
  labs(title="Dot plot", 
       subtitle="Drug treatment response",
       caption="Drug : mg",
       x="Drugs",
       y="Doses")+theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
######################################################################
### Pyramid bar plot 
options(scipen = 999)  # turns of scientific notations like 1e+40
# Read data
email_campaign_funnel <- read.csv("https://raw.githubusercontent.com/selva86/datasets/master/email_campaign_funnel.csv")
# X Axis Breaks and Labels 
brks <- seq(-15000000, 15000000, 5000000)
lbls = paste0(as.character(c(seq(15, 0, -5), seq(5, 15, 5))), "m")
# Plot
ggplot(email_campaign_funnel, aes(x = Stage, y = Users, fill = Gender)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(title="Email Campaign Funnel") +
  theme_tufte() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_brewer(palette = "Dark2")  # Color palette
