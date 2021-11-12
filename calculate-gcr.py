import pandas as pd

all_variants = []
all_genes = []
genes_with_gcr = []

variants_file = open('./27_Jun_condensed_ccg_with_vcr.csv', 'r')

for line in variants_file:
    if line.startswith("\ufeffGENE") == False:
        all_variants.append(line.strip("\n").split(","))

for variant in all_variants:
    if variant[0] not in all_genes:
        all_genes.append(variant[0])

for gene in all_genes:
    gene_gcr_inverse = 1
    gene_gcr_inverse_o = 1
    gene_gcr_inverse_la = 1
    gene_gcr_inverse_eur_fin = 1
    gene_gcr_inverse_am = 1
    gene_gcr_inverse_ea = 1
    gene_gcr_inverse_me = 1
    gene_gcr_inverse_sa = 1
    gene_gcr_inverse_aj = 1
    gene_gcr_inverse_aa = 1
    gene_gcr_inverse_eur_nonfin = 1
    for var in all_variants:
        if var[0] == gene:
            gene_gcr_inverse *= float(var[3])
            gene_gcr_inverse_o *= float(var[5])
            gene_gcr_inverse_la *= float(var[7])
            gene_gcr_inverse_eur_fin *= float(var[9])
            gene_gcr_inverse_am *= float(var[11])
            gene_gcr_inverse_ea *= float(var[13])
            gene_gcr_inverse_me *= float(var[15])
            gene_gcr_inverse_sa *= float(var[17])
            gene_gcr_inverse_aj *= float(var[19])
            gene_gcr_inverse_aa *= float(var[21])
            gene_gcr_inverse_eur_nonfin *= float(var[23])

    genes_with_gcr.append([gene, 1-gene_gcr_inverse_o, 1-gene_gcr_inverse_la, 1-gene_gcr_inverse_eur_fin, 1-gene_gcr_inverse_am, 1-gene_gcr_inverse_ea, 1-gene_gcr_inverse_me, 1-gene_gcr_inverse_sa, 1-gene_gcr_inverse_aj, 1-gene_gcr_inverse_aa, 1-gene_gcr_inverse_eur_nonfin, 1-gene_gcr_inverse])

df_genes_with_gcr = pd.DataFrame(genes_with_gcr, columns=['GENE', 'GCR OTHER', 'GCR LATINO/ADMIXED_AMERICAN', 'GCR EUROPEAN_(FINNISH)', 'GCR AMISH', 'GCR EAST_ASIAN', 'GCR MIDDLE_EASTERN', 'GCR SOUTH_ASIAN', 'GCR ASHKENAZI_JEWISH', 'GCR AFRICAN/AFRICAN-AMERICAN', 'GCR EUROPEAN_(NON-FINNISH)', 'GCR ALL'])

#print(df_genes_with_gcr)

df_genes_with_gcr.to_csv("28_Jun_Genes_with_GCR.csv", index=False)
