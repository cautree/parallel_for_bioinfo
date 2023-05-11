## density plots and histogram overlay
ggplot(covars, aes(x=hs_child_age_None)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#FF6666") +
 xlab("Child Age")
