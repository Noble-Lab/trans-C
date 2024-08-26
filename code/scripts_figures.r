
# --------------------------------------------------------------------------- #
# Fig 2 Plasmodium
# --------------------------------------------------------------------------- #

a <- rt("plasmodium")
ggplot(a) + geom_line(aes(x=V2, y=V1, color=V3), size=1) + 
xlab("False Positive Rate") + ylab("True Positive Rate") +
plot_helper() +
theme(legend.title = element_blank()) +
scale_colour_manual(values = c("#0072b2", "#E66100", "#5D3A9B"), labels=c("Random Walk", "Greedy")) +
theme(legend.text = element_text(color="black", size=22))

# 2B
a <- rt("plasmo_square")
ggplot(a) + geom_tile(aes(x=V1,y=V2,fill=V3)) +
plot_helper() + xlab("") + ylab("") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
theme(axis.ticks = element_blank()) + 
theme(axis.text = element_blank()) + 
theme(axis.title = element_text(color="black", size=16)) +
scale_fill_gradient2(low="#6B6B6B", mid="white",high="purple",midpoint = 0, name="reads",limits=c(-10, 30), oob=squish)

# 2C Baseline 
a <- rt("roc_plasmo_base")
ggplot(a) + geom_line(aes(x=V2, y=V1, color=V3, linetype=V4), size=1) + 
xlab("False Positive Rate") + ylab("True Positive Rate") +
plot_helper() +
theme(legend.title = element_blank()) +
scale_colour_manual(values = c("#5D3A9B", "#E66100"), labels=c("Schizont", "Trophozoite")) +
scale_linetype_manual(values = c("dotted", "solid"), labels=c("Greedy", "RW")) +
theme(legend.text = element_text(color="black", size=22))


# --------------------------------------------------------------------------- #
# Fig 3 Mouse Greek Islands  
# --------------------------------------------------------------------------- #
a <- rt("greek")
ggplot(a) + geom_line(aes(x=V2, y=V1, color=V3), size=1) + 
xlab("False Positive Rate") + ylab("True Positive Rate") +
plot_helper() +
theme(legend.title = element_blank()) +
scale_colour_manual(values = c("#0072b2", "#E66100", "#5D3A9B"), labels=c("mOSN", "HBC")) +
theme(legend.text = element_text(color="black", size=22))

# 3C
a <- rt("greek_average")
ggplot(a) + geom_tile(aes(x=V2,y=V1,fill=V3)) +
plot_helper() + xlab("") + ylab("") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
theme(axis.ticks = element_blank()) + 
theme(axis.text = element_blank()) + 
theme(axis.title = element_text(color="black", size=26)) +
scale_fill_gradient(low="white",high="purple",name="Average Weight",limits=c(0, 220), breaks=c(50,100,150,200))


# --------------------------------------------------------------------------- #
# Fig 4 RBM factory 
# --------------------------------------------------------------------------- #
a <- rt("roc_time_rbm")
ggplot(a) + geom_line(aes(x=V2, y=V1, color=V3), size=1) + 
xlab("False Positive Rate") + ylab("True Positive Rate") +
plot_helper() +
theme(legend.title = element_blank()) +
scale_colour_manual(values = c("#0072b2", "#E66100", "#5D3A9B"), labels=c("hPSCs", "Early CMs", "Late CMs")) +
theme(legend.text = element_text(color="black", size=22))

# 4D
a <- rt("rbm_square")
ggplot(a) + geom_tile(aes(x=V2,y=V1,fill=V3)) +
plot_helper() + xlab("") + ylab("") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
theme(axis.ticks = element_blank()) + 
theme(axis.text = element_blank()) + 
theme(axis.title = element_text(color="black", size=26)) +
scale_fill_gradient(low="white",high="purple",name="Average Weight",limits=c(0, 320), breaks=c(10,20))


# --------------------------------------------------------------------------- #
# Fig 5 TFs
# --------------------------------------------------------------------------- #
# 5B
a <- rt("weights")
ggplot(a) + geom_violin(aes(x=V2, y=V1, fill=V2)) +
xlab("") + ylab("subnetwork weight") +
theme_bobby_white() +  scale_y_continuous(trans='log2') + 
scale_x_discrete(labels = c("random\nbins", "random\nseed", "TF loci", "matched\nrandom\nseed", "TF-based\nseed")) +
scale_fill_manual(labels = c("random bins", "random seed", "TF loci", "matched random seed", "TF-based seed"), values= c("#009e73", "#d5c711", "#0072b2", "#E66100", "#5D3A9B")) +
theme(legend.position="none") 

# 5C
a <- rt("dbp")
ggplot() +
plot_helper() + xlab("matched random subnetwork weight") + ylab("DBP-based subnetwork weight") +
geom_errorbar(data=eb, aes(y = V2, xmin = V3, xmax = V4, width=0.002), color="grey") +
geom_point(data=a, aes(y=V2,x=V3, color=V4), size=rel(2.3)) +
geom_abline(intercept = 0, color="black") +
scale_y_continuous(trans='log2') + 
scale_x_continuous(trans='log2') +
scale_color_manual(values = c("fdr"="#c3352b", "no"="#352bc3")) +
geom_text(data=a, aes(y=V2,x=V3,label=V1), hjust=0, vjust=0) +
theme(legend.position="none")

# 5 E
a <- rt("tf_son")
ggplot(a) + geom_dot(aes(x=V2,y=V1,fill=V3)) +
xlab("Top ranked loci (DBP clique)") + ylab("Bottom ranked loci") +
plot_helper() +
geom_abline(intercept = 0, color="black", linetype="dotted")

# 5 G
a <- rt("tf_son_sig")
ggplot(a) + geom_violin(aes(x=V2,y=V1,fill=V2)) +
xlab("") + ylab("DBP intrinsic disorder score") +
plot_helper() + theme(legend.position="none") +
coord_flip()

# 5 H
a <- rt("idp_scores")
ggplot(a) + geom_violin(aes(x=V2,y=V1,fill=V2)) +
xlab("") + ylab("DBP intrinsic disorder score") +
plot_helper() + theme(legend.position="none") +
coord_flip()


# --------------------------------------------------------------------------- #
# Fig 6 
# --------------------------------------------------------------------------- #
a <- rt("eclip")
ggplot() + geom_point(data=a, aes(x=V2,y=V3, color=V4)) +
plot_helper() + xlab("RBP-based subnetwork") + ylab("Matched random subnetwork") +
geom_errorbar(data=eb, aes(x = V2, ymin = V3, ymax = V4, width=0.02, color=V5)) +
geom_abline(intercept = 0, color="black") +
scale_y_continuous(trans='log2') + 
scale_x_continuous(trans='log2') +
scale_color_manual(values = c("fdr"="#c3352b", "no"="#352bc3")) +
theme(legend.position="none")
# 6B
a <- rt("tf_son")
ggplot(a) + geom_dot(aes(x=V2,y=V1,fill=V3)) +
xlab("Top ranked loci (RBP clique)") + ylab("Bottom ranked loci") +
plot_helper() +
geom_abline(intercept = 0, color="black", linetype="dotted")

# --------------------------------------------------------------------------- #
# Fig S2 
# --------------------------------------------------------------------------- #
a <- rt("roc_high_conf")
ggplot(a) + geom_line(aes(x=V2, y=V1, color=V3), size=1) + 
xlab("False Positive Rate") + ylab("True Positive Rate") +
plot_helper() +
theme(legend.title = element_blank()) +
scale_colour_manual(values = c("#5D3A9B"), labels=c("high\nconfidence")) +
theme(legend.text = element_text(color="black", size=22))


# --------------------------------------------------------------------------- #
# Fig S3 left ventricle
# --------------------------------------------------------------------------- #
a <- rt("left_ventricle")
ggplot(a) + geom_tile(aes(x=V2,y=V1,fill=V3)) +
plot_helper() + xlab("") + ylab("") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
theme(axis.ticks = element_blank()) + 
theme(axis.text = element_blank()) + 
theme(axis.title = element_text(color="black", size=26)) +
scale_fill_gradient(low="white",high="purple",name="Pearson Correlation")


# --------------------------------------------------------------------------- #
# Fig S4
# --------------------------------------------------------------------------- #
a <- rt("hypergeometric")
ggplot(a) + geom_tile(aes(x=V2,y=V1,fill=V3)) +
plot_helper() + xlab("DPB-seeded subnetwork") + ylab("DNA binding protein") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
theme(axis.ticks = element_blank()) + 
theme(axis.text = element_blank()) + 
theme(axis.title = element_text(color="black", size=26)) +
scale_fill_gradient(low="white",high="purple",name="-log(p-value)",limits=c(0, 70), breaks=c(20,40,60), oob=squish)

# S4B
a <- rt("permuted")
ggplot() + geom_point(data=a, aes(x=V2,y=V3)) +
plot_helper() + xlab("Matched random subnetwork weight\nin a randomly permuted network") + ylab("DBP-based subnetwork weight") +
geom_errorbar(data=eb, aes(x = V2, ymin = V3, ymax = V4, width=0.02)) +
geom_abline(intercept = 0, color="black") +
scale_y_continuous(trans='log2') + 
scale_x_continuous(trans='log2') +
theme(legend.position="none")

# S4C
a <- rt("top_v_bottom")
ggplot(a) + geom_jitter(aes(x=V3, y=V1, color=V3)) + 
theme_bw() + xlab("") + ylab("preference to bind A compartment") +
theme(axis.text = element_text(color="black", size=14)) + 
theme(axis.title = element_text(color="black", size=15)) +
theme(legend.text=element_text(color="black", size=15)) +
theme(legend.title = element_blank()) +
scale_colour_manual(values = c("#5D3A9B", "#E66100")) +
theme(legend.position="none")


# --------------------------------------------------------------------------- #
# S5 SPRITE
# --------------------------------------------------------------------------- #
a <- rt("sprite_weight")
ggplot(a) + geom_violin(aes(x=V3, y=V2, fill=V3)) +
xlab("") + ylab("proximity") +
theme_bobby_white() +  
scale_x_discrete(labels = c("random\nbins", "trans-C\nclique")) +
scale_fill_manual(labels = c("random\nbins", "trans-C\nclique"), values= c("#009e73", "#d5c711")) +
theme(legend.position="none")
# S5B
a <- rt("sprite")
ggplot() +
plot_helper() + xlab("matched random SPRITE subnetwork weight") + ylab("DBP-based subnetwork weight") +
geom_errorbar(data=eb, aes(y = V2, xmin = V3, xmax = V4, width=0.002), color="grey") +
geom_point(data=a, aes(y=V2,x=V3, color=V4), size=rel(2.3)) +
geom_abline(intercept = 0, color="black") +
scale_y_continuous(trans='log2') + 
scale_x_continuous(trans='log2') +
scale_color_manual(values = c("fdr"="#c3352b", "no"="#352bc3")) +
theme(legend.position="none")
a <- rt("sprite_fold")
ggplot(a) + geom_point(aes(x=V2, y=V3, color=V4)) +
plot_helper() + xlab("fold change over matched random (Hi-C)") + ylab("fold change over matched random (SPRITE)") +
geom_abline(intercept = 0, color="black") +
scale_color_manual(values = c("#cd00cd", "#008000", "#ffa500", "#000000"), labels = c("both", "Hi-C only", "Sprite only", "neither"), name="significant") +
scale_y_continuous(breaks=c(0.5,1,1.5,2.0), labels = c("0.5","1","1.5","2.0")) +
coord_cartesian(ylim=c(0.5,2)) +
theme(legend.position = c(.82,.85))


# --------------------------------------------------------------------------- #
# S6 MERFISH
# --------------------------------------------------------------------------- #
a <- rt("dist")
ggplot(a) + geom_violin(aes(x=V3, y=V2, fill=V3)) +
xlab("") + ylab("proximity") +
theme_bobby_white() +  
scale_x_discrete(labels = c("random\nbins", "trans-C\nclique")) +
scale_fill_manual(labels = c("random\nbins", "trans-C\nclique"), values= c("#009e73", "#d5c711")) +
theme(legend.position="none")
a <- rt("bottom")
ggplot(a) + geom_point(aes(x=V2, y=V1)) +
xlab("proximity of trans-C bottom ranked loci") + ylab("proximity of trans-C top ranked loci") +
geom_abline(intercept = 0, color="black") +
coord_cartesian(xlim = c(0.005, 0.012), ylim = c(0.005, 0.012)) +
theme_bobby_white() 
ggplot() +
geom_errorbar(data=eb, aes(y = V2, xmin = V3, xmax = V4, width=0.002), color="grey") +
geom_point(data=a, aes(y=V2,x=V3, color=V4), size=rel(2.3)) +
geom_abline(intercept = 0, color="black") +
scale_color_manual(values = c("fdr"="#c3352b", "no"="#352bc3")) +
theme(legend.position="none")


# --------------------------------------------------------------------------- #
# helper function
# --------------------------------------------------------------------------- #

plot_helper <- function (base_size = 12, base_family = "") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace% 
        theme(
            axis.text = element_text(color="black", size=16),
            axis.title = element_text(color="black", size=17),
            legend.text = element_text(color="black", size=18),
			legend.title = element_text(color="black", size=18),
			panel.grid.minor = element_blank()
    )   
}
