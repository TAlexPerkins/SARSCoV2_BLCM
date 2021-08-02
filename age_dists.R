## Load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(data.table,
       viridis)

## Read in demographic data
ages <- fread("NDID_Demographics.csv")

## Plotting controls
inset <- TRUE ## inset plot of upper values
groups <- FALSE ## color students & staff differently

## Plot histogram of ages
pdf(ifelse(groups,
    ifelse(inset,"age_distribution.pdf","age_distribution_no_inset.pdf"),
    ifelse(inset,"age_distribution_no_groups.pdf","age_distribution_no_groups_or_inset.pdf")),
    height=ifelse(inset,4,3.25),width=ifelse(inset,6,3.25))
par(mar=c(4.1,4.1,0.6,0.6))
hist(ages$AGE,breaks=seq(min(ages$AGE)-0.5,max(ages$AGE)+0.5),
     col=ifelse(groups,viridis(10)[1],"lightgray"),yaxs="i",xlab="Age",
     ylab="Number of participants",main="")
if (groups) {
    hist(ages$AGE[ages$DESCRIPTION_1=="Student"],add=TRUE,
         breaks=seq(min(ages$AGE)-0.5,max(ages$AGE)+0.5),
         col=viridis(10)[6],yaxs="i")
}
if (inset) {
    par(fig=c(0.3,0.9,0.2,0.9),new=TRUE)
    min.age <- ifelse(groups,min(ages$AGE[ages$DESCRIPTION_1=="Staff"]),30)
    hist(ages$AGE[ages$AGE>=min.age],breaks=seq(min.age-0.5,max(ages$AGE)+0.5),
         col=ifelse(groups,viridis(10)[1],"lightgray"),yaxs="i",xlab="",ylab="",main="")
    if (groups) {
        hist(ages$AGE[ages$AGE>=min.age & ages$DESCRIPTION_1=="Student"],add=TRUE,
             breaks=seq(min.age-0.5,max(ages$AGE)+0.5),
             col=viridis(10)[6],yaxs="i")
        legend("right",bty="n",fill=viridis(10)[c(6,1)],legend=c("Student","Staff"))
    }
}
dev.off()

## Print summary statistics to screen
summary(ages$AGE)
summary(ages$AGE[ages$DESCRIPTION_1=="Student"])
summary(ages$AGE[ages$DESCRIPTION_1=="Staff"])
## Percentage students between 18 and 22
sum(ages$AGE[ages$DESCRIPTION_1=="Student"]  <= 22)/nrow(ages)
## Percentage staff 65 or older
sum(ages$AGE[ages$DESCRIPTION_1=="Staff"]  >= 65)
