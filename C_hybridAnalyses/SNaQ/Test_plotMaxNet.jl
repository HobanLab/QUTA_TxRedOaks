# Specify filepath to the starting topology (a Newick text file), and read in
SNaQ_AllTips_15sp_MaxNet_filepath = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/AllTips1/Output/Net1_AllTips_15sp_MaxNet.tre"
SNaQ_AllTips_MaxNet = readTopology(SNaQ_AllTips_15sp_MaxNet_filepath)

# Plot the resulting phylogenetic network and save to a PDF file
using RCall
R"pdf"("Plot-AllTips_Net1.pdf", width=15, height=9);
plot(SNaQ_AllTips_MaxNet, :R);
R"dev.off()";
