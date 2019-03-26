library(ggplot2)
library(dplyr)

sim_df=readRDS("sim_df.rds")

# ORCHARD MAP 
orchard1%>%filter(!subsppcyt=="AR")%>%
  ggplot(aes(x=x,y=y,size=(volume_2017)*10^-6,colour=subsppcyt))+geom_point()+
  #geom_jitter(aes(y=size_t,x=subsppcyt, alpha=0.35,colour="#05C3DEBF"))+
  xlab("x")+ylab("y")+
  ggtitle(expression("Orchard common garden:"~italic(A.~tridentata)))+theme_bw()+
  #scale_fill_manual(name = "Crown volume", guide = "legend",labels = "") +
  #scale_colour_manual(name = 'Suspecies:cytotype',values=gg_color_hue(5), labels = c('T2n','T4n','V2n','V4n','W4n'))
  guides(size=guide_legend("Crown volume",labels=""),colour=guide_legend("Subspecies:cytotype"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust=0))


#
ggplot(df,aes(x=spacing,y=rel_cov))+geom_point()+
  theme_bw()+labs(x="Spacing (m)",y="Relative cover (%)")+ggtitle("")


ggplot(df,aes(x=spacing,y=palive))+geom_point()+
  theme_bw()+labs(x="Spacing (m)",y="Survival (%)")+ggtitle("")
  

ggplot(df)+geom_point(aes(x=spacing,y=avg_height),colour="#C77CFF")+
  theme_bw()+labs(x="Spacing (m)",y="Average height (m)")+ggtitle("")+
  theme(legend.position ="none",panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust=0))
 

df%>%filter(rel_cov<1)%>%
ggplot()+
  geom_point(aes(y=rel_cov,x=spacing,colour="#F8766D"),shape=19)+
  geom_point(aes(y=palive,x=spacing,colour="#7CAE00"),shape=19)+
  theme_bw()+labs(x="Spacing (m)",y="Proportion (%)")+ggtitle("")+
  scale_colour_manual(name="", values=c("#F8766D","#7CAE00"),labels=c("Canopy cover","Survival"))+
  theme(legend.position="top",panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust=0))+
  guides(colour = guide_legend(override.aes = list(shape = c(19,19))))
  
plot(df$rel_cov~df$spacing, ylab="", xlab="Spacing (m)",col="purple",pch=19)
plot(df$palive~df$spacing, col="red",pch=19)
plot(df$avg_height~df$spacing, col="blue",pch=19)
#legend("right",legend=c("# alive individuals","% cover"),col=c("red","blue"),pch=1,lty=2,lwd=3,bty="n")
