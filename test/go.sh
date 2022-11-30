bluues2 protein.pqr protein_out  
bluues2 protein.pqr protein_out -srf
bluues2 protein.pqr protein_out -srfpot
bluues2 protein.pqr protein_out -dx
bluues2 protein.pqr protein_out -pka
bluues2 dna.pqr dna_out -pka -pkadef pkadefs.txt 
vmd -e gbr.vmd
vmd -e heatmap.vmd
vmd -e isosurface.vmd
vmd -e srfpot.vmd
