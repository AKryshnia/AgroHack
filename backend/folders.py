import os

# Определяем путь к директории, где находится скрипт
script_dir = os.path.dirname(os.path.abspath(__file__))

data_dir = os.path.join(script_dir, "data")
output_dir = os.path.join(script_dir, "output")

# Пути к файлам
ancient_barley_nuc_file = os.path.join(data_dir, "ancient_barley.fasta")
vulgare_file = os.path.join(data_dir, "hordeum_vulgare.fasta")
spontaneum_file = os.path.join(data_dir, "spontaneum_genes.fasta")
reference_genes = ["Btr1", "Btr2", "Vrs1"]
reference_files = [
    os.path.join(data_dir, f"{gene}_reference.fasta") for gene in reference_genes
]

combined_nuc_file = os.path.join(output_dir, "combined_nucleotide_sequences.fasta")
aligned_nuc_file = os.path.join(output_dir, "aligned_nucleotide_sequences.fasta")
combined_prot_file = os.path.join(output_dir, "combined_protein_sequences.fasta")
aligned_prot_file = os.path.join(output_dir, "aligned_protein_sequences.fasta")

differences_file = os.path.join(output_dir, "differences_table.txt")
similarity_file = os.path.join(output_dir, "similarities_table.txt")
consensus_file = os.path.join(output_dir, "consensus_sequences.txt")
amino_acid_frequency_file = os.path.join(output_dir, "amino_acid_frequencies.png")
mutation_analysis_file = os.path.join(output_dir, "mutation_analysis.xlsx")
morphology_report_file = os.path.join(output_dir, "morphology_report.txt")
phylogenetic_tree_file = os.path.join(output_dir, "phylogenetic_tree.png")
region_report_file = os.path.join(output_dir, "region_report.txt")
