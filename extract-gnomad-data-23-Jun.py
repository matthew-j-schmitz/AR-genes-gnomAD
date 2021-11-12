import pandas as pd
import numpy as np
import os

# will take duplicate of clinvar variant file, create array,

combined_clinvar_gnomad_list = []

datapath = "./gnomad-gene-csv-files"

def generate_clinvar_list(clinvar_file):
    # Generates a list of lists of clinvar variants from csv file

    # Open csv file of plp clinvar variants
    clinvar_f = open(clinvar_file, 'r')

    # Generate empty list to store clinvar variants
    list_of_clinvar_variants = []

    # Add each line in the csv as an entry in the list
    for line in clinvar_f:
        parsed_line = line.strip("\n").split("\t")
        parsed_line.pop(0)
        list_of_clinvar_variants.append(parsed_line)

    # Remove the header line from the list
    list_of_clinvar_variants.pop(0)

    # Close csv file of plp clinvar variants
    clinvar_f.close()

    return list_of_clinvar_variants

def generate_gnomad_list(gnomad_file):
    # Generates a list of lists of clinvar variants from gnomad csv file

    gnomad_f = open(gnomad_file, 'r')

    list_of_gnomad_variants = []

    # Add each line in the csv as an entry in the list
    for line in gnomad_f:
        parsed_line = line.strip("\n").split(",")
        list_of_gnomad_variants.append(parsed_line)

    # Remove the header line from the list
    list_of_gnomad_variants.pop(0)

    gnomad_f.close()

    return list_of_gnomad_variants

def extract_gnomad_data_by_variant_ID(gnomad_list, clinvar_list, combined_list):
    # Extracts data from gnomad list and adds it to clinvar list

    # To store index of located variants
    found_variant_indices = []

    # Iterate through clinvar list and identify presence in gnomad list
    for i in range(len(clinvar_list)):
        for gnomad_variant in gnomad_list:
            if clinvar_list[i][-1] == gnomad_variant[0] + "-" + gnomad_variant[1] + "-" + gnomad_variant[3] + "-" + gnomad_variant[4]:
                clinvar_line = clinvar_list[i]
                for n in list(range(-45, 0)):
                    clinvar_line.append(gnomad_variant[n])
                combined_list.append(clinvar_line)
                found_variant_indices.append(i)

    return found_variant_indices

my_clinvar_variants = generate_clinvar_list("./10_Jun_selected_autorec_plp_nonconf_tsv.tsv")

#print(len(my_clinvar_variants))

gene_count = 0

# Iterate through gnomad csv files
for single_gene_file in os.listdir(datapath):

    gene_count += 1

    if single_gene_file.startswith("gnomAD"):

        file_location = "./gnomad-gene-csv-files/" + single_gene_file

        current_file_gnomad_list = generate_gnomad_list(file_location)

        indices_of_found_variants = extract_gnomad_data_by_variant_ID(current_file_gnomad_list, my_clinvar_variants, combined_clinvar_gnomad_list)

        # Remove the variant from the clinvar list after it has been found
        for index in sorted(indices_of_found_variants, reverse=True):
            my_clinvar_variants.pop(index)

    # Keep track of how many genes have been completed
    print(gene_count)

print(len(combined_clinvar_gnomad_list))
print(len(my_clinvar_variants))

#print(combined_clinvar_gnomad_list[0:3])

gnomad_clinvar_df = pd.DataFrame(combined_clinvar_gnomad_list, columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'ALLELEID', 'CLNDISDB', 'CLNHGVS', 'CLNSIG', 'CLNVC', 'MC', 'AF_ESP', 'AF_EXAC', 'AF_TGP', 'RS', 'CLNSIGCONF', 'CLNDISDBINCL', 'GENEINFO', 'GENE', 'GNOMAD_ID', "ALLELE_COUNT", "ALLELE_NUMBER", "ALLELE_FREQUENCY", "HOMOZYGOTE_COUNT", "HEMIZYGOTE_COUNT", "ALLELE_COUNT_OTHER", "ALLELE_NUMBER_OTHER", "HOMOZYGOTE_COUNT_OTHER", "HEMIZYGOTE_COUNT_OTHER", "ALLELE_COUNT_LATINO/ADMIXED_AMERICAN", "ALLELE_NUMBER_LATINO/ADMIXED_AMERICAN", "HOMOZYGOTE_COUNT_LATINO/ADMIXED_AMERICAN", "HEMIZYGOTE_COUNT_LATINO/ADMIXED_AMERICAN", "ALLELE_COUNT_EUROPEAN_(FINNISH)", "ALLELE_NUMBER_EUROPEAN_(FINNISH)", "HOMOZYGOTE_COUNT_EUROPEAN_(FINNISH)", "HEMIZYGOTE_COUNT_EUROPEAN_(FINNISH)", "ALLELE_COUNT_AMISH", "ALLELE_NUMBER_AMISH", "HOMOZYGOTE_COUNT_AMISH", "HEMIZYGOTE_COUNT_AMISH", "ALLELE_COUNT_EAST_ASIAN", "ALLELE_NUMBER_EAST_ASIAN", "HOMOZYGOTE_COUNT_EAST_ASIAN", "HEMIZYGOTE_COUNT_EAST_ASIAN", "ALLELE_COUNT_MIDDLE_EASTERN", "ALLELE_NUMBER_MIDDLE_EASTERN", "HOMOZYGOTE_COUNT_MIDDLE_EASTERN", "HEMIZYGOTE_COUNT_MIDDLE_EASTERN", "ALLELE_COUNT_SOUTH_ASIAN", "ALLELE_NUMBER_SOUTH_ASIAN", "HOMOZYGOTE_COUNT_SOUTH_ASIAN", "HEMIZYGOTE_COUNT_SOUTH_ASIAN", "ALLELE_COUNT_ASHKENAZI_JEWISH", "ALLELE_NUMBER_ASHKENAZI_JEWISH", "HOMOZYGOTE_COUNT_ASHKENAZI_JEWISH", "HEMIZYGOTE_COUNT_ASHKENAZI_JEWISH", "ALLELE_COUNT_AFRICAN/AFRICAN-AMERICAN", "ALLELE_NUMBER_AFRICAN/AFRICAN-AMERICAN", "HOMOZYGOTE_COUNT_AFRICAN/AFRICAN-AMERICAN", "HEMIZYGOTE_COUNT_AFRICAN/AFRICAN-AMERICAN", "ALLELE_COUNT_EUROPEAN_(NON-FINNISH)", "ALLELE_NUMBER_EUROPEAN_(NON-FINNISH)", "HOMOZYGOTE_COUNT_EUROPEAN_(NON-FINNISH)", "HEMIZYGOTE_COUNT_EUROPEAN_(NON-FINNISH)"])

print("gnomAD df:")
print(gnomad_clinvar_df)

clinvar_remaining_df = pd.DataFrame(my_clinvar_variants, columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'ALLELEID', 'CLNDISDB', 'CLNHGVS', 'CLNSIG', 'CLNVC', 'MC', 'AF_ESP', 'AF_EXAC', 'AF_TGP', 'RS', 'CLNSIGCONF', 'CLNDISDBINCL', 'GENEINFO', 'GENE', 'GNOMAD_ID'])

print("Clinvar remaining df:")
print(clinvar_remaining_df)

gnomad_clinvar_df.to_csv("23_Jun_combined_clinvar_gnomad.csv")

gnomad_clinvar_df.to_csv("23_Jun_combined_clinvar_gnomad_tsv.tsv", sep="\t")

clinvar_remaining_df.to_csv("23_Jun_remaining_clinvar.csv")

clinvar_remaining_df.to_csv("23_Jun_remaining_clinvar.tsv", sep="\t")

# Took 1 min 38 s for 5 gnomad files (97 variants)
# Took 195 s for 13 gnomad files (294 variants), 3 min 11 s with removal
# Took 405 s for 30 gnomad files (610 variants) with removal - 0.66 s / var
# Took 700 s for 57 gnomad files (1326 variants) with removal - 0.53 s / var
# Took approximately 1.5 hours for all 419 gnomad files (8437 variants)
