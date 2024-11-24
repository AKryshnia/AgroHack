# barley_analysis.py

import os
import logging
from Bio import SeqIO, Entrez, AlignIO
import matplotlib.pyplot as plt

from utils import (
    create_dir_if_not_exists,
    check_clustal_omega_installed,
    calculate_amino_acid_frequency,
    export_mutations_to_excel,
    download_sequence_by_accession,
    download_sequences_for_genes,
    filter_sequences_by_genes,
    run_clustal_omega,
    identify_advanced_mutations,
    process_sequences,
    identify_amino_acid_changes,
    get_consensus_sequence,
    reconstruct_morphology,
    compute_differences_and_similarity,
    perform_phylogenetic_analysis,
)
from folders import (
    data_dir,
    output_dir,
    script_dir,
    ancient_barley_nuc_file,
    vulgare_file,
    spontaneum_file,
    combined_nuc_file,
    aligned_nuc_file,
    combined_prot_file,
    aligned_prot_file,
    differences_file,
    similarity_file,
    consensus_file,
    amino_acid_frequency_file,
)

# Укажите свою почту для работы с NCBI
Entrez.email = "your_email@example.com"

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s: %(message)s",
    filename="analysis.log",
)
logger = logging.getLogger(__name__)


create_dir_if_not_exists(data_dir)
create_dir_if_not_exists(output_dir)
create_dir_if_not_exists(script_dir)


# Главная функция
def main():
    # Проверяем наличие директории data
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Проверяем установку Clustal Omega
    check_clustal_omega_installed()

    genes = ["Btr1", "Btr2", "Vrs1"]

    # Скачиваем данные, если файлы не существуют
    # Скачивание референсных последовательностей
    ref_accession_numbers = {"Btr1": "AK373702", "Btr2": "AK368246", "Vrs1": "AY093351"}
    reference_sequences = []
    for gene, acc_num in ref_accession_numbers.items():
        output_file = os.path.join(data_dir, f"{gene}_reference.fasta")
        if not os.path.exists(output_file):
            success = download_sequence_by_accession(acc_num, output_file)
            if not success:
                logger.error(
                    f"Не удалось скачать референсную последовательность {gene}"
                )
                continue

        record = SeqIO.read(output_file, "fasta")
        record.id = gene + "_ref"
        reference_sequences.append(record)

    if not reference_sequences:
        logger.error(
            "Нет доступных референсных последовательностей. Прерывание анализа."
        )
        return

    # Ограничиваем количество и размер скачиваемых последовательностей
    max_results = 5
    max_sequence_length = 2000

    if not os.path.exists(vulgare_file):
        download_sequences_for_genes(
            genes,
            "Hordeum vulgare",
            vulgare_file,
            max_results=max_results,
            max_sequence_length=max_sequence_length,
        )

    if not os.path.exists(spontaneum_file):
        download_sequences_for_genes(
            genes,
            "Hordeum vulgare subsp. spontaneum",
            spontaneum_file,
            max_results=max_results,
            max_sequence_length=max_sequence_length,
        )

    # Чтение древних последовательностей
    if not os.path.exists(ancient_barley_nuc_file):
        logger.error(
            f"Файл {ancient_barley_nuc_file} не найден. Пожалуйста, поместите его в директорию data."
        )
        return

    ancient_nuc_sequences = list(SeqIO.parse(ancient_barley_nuc_file, "fasta"))
    ancient_nuc_sequences = filter_sequences_by_genes(ancient_nuc_sequences, genes)
    for seq in ancient_nuc_sequences:
        seq.id = "Ancient_" + seq.id

    # Чтение современных последовательностей
    vulgare_nuc_sequences = list(SeqIO.parse(vulgare_file, "fasta"))
    vulgare_nuc_sequences = filter_sequences_by_genes(vulgare_nuc_sequences, genes)
    for seq in vulgare_nuc_sequences:
        seq.id = "Vulgare_" + seq.id

    spontaneum_nuc_sequences = list(SeqIO.parse(spontaneum_file, "fasta"))
    spontaneum_nuc_sequences = filter_sequences_by_genes(
        spontaneum_nuc_sequences, genes
    )
    for seq in spontaneum_nuc_sequences:
        seq.id = "Spontaneum_" + seq.id

    # Объединение всех нуклеотидных последовательностей
    all_nuc_sequences = (
        ancient_nuc_sequences
        + vulgare_nuc_sequences
        + spontaneum_nuc_sequences
        + reference_sequences
    )
    if not all_nuc_sequences:
        logger.error("Нет последовательностей для анализа. Прерывание программы.")
        return

    SeqIO.write(all_nuc_sequences, combined_nuc_file, "fasta")
    logger.info(
        f"Объединенные нуклеотидные последовательности сохранены в {combined_nuc_file}"
    )

    # Выравнивание нуклеотидных последовательностей (Критерий 2)
    run_clustal_omega(combined_nuc_file, aligned_nuc_file, seqtype="DNA")
    alignment_nuc = AlignIO.read(aligned_nuc_file, "fasta")
    logger.info(
        f"Количество выравненных нуклеотидных последовательностей: {len(alignment_nuc)}"
    )

    # Идентификация мутаций (Критерий 1)
    mutations, mutation_stats = identify_advanced_mutations(
        all_nuc_sequences, reference_sequences[0]
    )
    export_mutations_to_excel(mutations)

    # Чтение и обработка белковых последовательностей
    # Трансляция нуклеотидных последовательностей в белковые
    all_protein_sequences = process_sequences(all_nuc_sequences)
    SeqIO.write(all_protein_sequences, combined_prot_file, "fasta")
    logger.info(
        f"Объединенные белковые последовательности сохранены в {combined_prot_file}"
    )

    # Выравнивание белковых последовательностей (Критерий 3)
    run_clustal_omega(combined_prot_file, aligned_prot_file, seqtype="Protein")
    alignment_prot = AlignIO.read(aligned_prot_file, "fasta")
    logger.info(
        f"Количество выравненных белковых последовательностей: {len(alignment_prot)}"
    )

    # Анализ замен аминокислот
    amino_acid_changes = identify_amino_acid_changes(
        alignment_prot, reference_sequences[0].id
    )
    with open(
        os.path.join(output_dir, "amino_acid_changes.txt"), "w", encoding="utf-8"
    ) as f:
        for seq_id, changes in amino_acid_changes.items():
            f.write(
                f"Замены аминокислот в последовательности {seq_id} относительно {reference_sequences[0].id}:\n"
            )
            for pos, ref_aa, aa in changes:
                f.write(f"Позиция {pos}: {ref_aa} -> {aa}\n")
            f.write("\n")
    logger.info("Замены аминокислот сохранены в файл: amino_acid_changes.txt")

    # Получение консенсусной последовательности
    consensus, degenerate_consensus = get_consensus_sequence(alignment_prot)
    with open(consensus_file, "w", encoding="utf-8") as f:
        f.write("Консенсусная последовательность:\n")
        f.write(str(consensus) + "\n\n")
        f.write("Дегенеративная консенсусная последовательность:\n")
        f.write(str(degenerate_consensus) + "\n")
    logger.info(
        "Консенсусные последовательности сохранены в файл: consensus_sequences.txt"
    )

    # Частота аминокислот
    frequencies = calculate_amino_acid_frequency(alignment_prot)

    # Построение гистограммы
    amino_acids = sorted(frequencies.keys())
    freq_values = [frequencies[aa] for aa in amino_acids]
    plt.figure(figsize=(12, 6))
    plt.bar(amino_acids, freq_values, color="skyblue")
    plt.xlabel("Аминокислоты")
    plt.ylabel("Частота")
    plt.title("Частота аминокислот в выравненных последовательностях")
    plt.tight_layout()
    plt.savefig(amino_acid_frequency_file)
    logger.info("График частоты аминокислот сохранен.")

    # Подсчет различий и схожести между группами (Критерий 5)
    # Разделение последовательностей по группам
    ancient_seqs = [
        seq for seq in all_protein_sequences if seq.id.startswith("Ancient_")
    ]
    vulgare_seqs = [
        seq for seq in all_protein_sequences if seq.id.startswith("Vulgare_")
    ]
    spontaneum_seqs = [
        seq for seq in all_protein_sequences if seq.id.startswith("Spontaneum_")
    ]

    # Подсчет различий и схожести
    differences, similarity = compute_differences_and_similarity(
        ancient_seqs, vulgare_seqs, spontaneum_seqs
    )

    # Сохранение результатов
    with open(differences_file, "w", encoding="utf-8") as f:
        f.write("Подсчет различий между группами:\n")
        for key, value in differences.items():
            f.write(f"{key}: {value} замен\n")
    logger.info("Таблица с подсчетом различий сохранена.")

    with open(similarity_file, "w", encoding="utf-8") as f:
        f.write("Схожесть между группами:\n")
        for key, value in similarity.items():
            f.write(f"{key}: {value:.2f}% схожести\n")
    logger.info("Таблица с подсчетом схожести сохранена.")

    # Реконструкция морфологии (Критерий 4)
    reconstruct_morphology(mutations)

    # Филогенетический анализ (Критерий 5)
    perform_phylogenetic_analysis(alignment_prot)

    logger.info("Анализ завершен.")


if __name__ == "__main__":
    main()
