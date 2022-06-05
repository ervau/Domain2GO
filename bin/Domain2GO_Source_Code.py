import pandas as pd

# Loading InterPro annotation data

protein2ipr = pd.read_csv('unique_ids.dat', names=["Uniprot", "Accession"], delimiter=' ')
ipr_domains = pd.read_csv('interpro_domains.tsv', sep='\t')
protein2domain = protein2ipr.merge(ipr_domains, on=['Accession'])
protein2domain.to_csv('domain_uniprot.txt', sep=' ', index=False)

# Loading GO annotation data with manual evidence code
go_annot_manual = pd.read_csv('manual_GOA_20191024_propagated.tsv', delimiter = "\t")
go_annot_manual = go_annot_manual[["GO_ID", "DB_OBJECT_ID"]]
go_annot_manual_unique=go_annot_manual.drop_duplicates()
go_annot_manual_unique.to_csv('goa_unique.txt', sep=' ', index=False)

# 