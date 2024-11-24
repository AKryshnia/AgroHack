# barley_analysis.py

import os
import sys
import shutil
import logging
import subprocess
import time
from collections import Counter
from Bio import SeqIO, Entrez, AlignIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import pandas as pd

# Укажите свою почту для работы с NCBI
Entrez.email = "your_email@example.com"

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s: %(message)s',
    filename='analysis.log'
)
logger = logging.getLogger(__name__)

# Определяем путь к директории, где находится скрипт
script_dir = os.path.dirname(os.path.abspath(__file__))

data_dir = os.path.join(script_dir, "data")
output_dir = os.path.join(script_dir, "output")

# Пути к файлам
ancient_barley_nuc_file = os.path.join(data_dir, "ancient_barley.fasta")
vulgare_file = os.path.join(data_dir, "hordeum_vulgare.fasta")
spontaneum_file = os.path.join(data_dir, "spontaneum_genes.fasta")
reference_genes = ["Btr1", "Btr2", "Vrs1"]
reference_files = [os.path.join(data_dir, f"{gene}_reference.fasta") for gene in reference_genes]

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


def create_dir_if_not_exists(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
        logger.info(f"Директория {dir} была создана.")
    else:
        logger.info(f"Директория {dir} уже существует.")


create_dir_if_not_exists(data_dir)
create_dir_if_not_exists(output_dir)
create_dir_if_not_exists(script_dir)


# Функция для проверки установки Clustal Omega
def check_clustal_omega_installed():
    if not shutil.which("clustalo"):
        raise EnvironmentError(
            "Clustal Omega не установлен. Убедитесь, что он доступен в PATH."
        )


# Функция для выполнения Clustal Omega
def run_clustal_omega(input_file, output_file, seqtype="Protein"):
    try:
        subprocess.run(
            [
                "clustalo",
                "-i",
                input_file,
                "-o",
                output_file,
                "--auto",
                f"--seqtype={seqtype}",
                "--force",
            ],
            check=True,
        )
        logger.info(f"Выравнивание завершено. Результаты сохранены в {output_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Ошибка при выполнении Clustal Omega: {e}")
        sys.exit(1)


# Функция для подсчета частоты аминокислот
def calculate_amino_acid_frequency(alignment):
    amino_acid_counts = Counter()
    for record in alignment:
        amino_acid_counts.update(
            str(record.seq).replace("-", "")
        )  # Убираем gap-символы
    total_count = sum(amino_acid_counts.values())
    frequencies = {aa: count / total_count for aa, count in amino_acid_counts.items()}
    return frequencies


# Функция для скачивания последовательности по Accession Number
def download_sequence_by_accession(accession_number, output_file, max_retries=3):
    for attempt in range(max_retries):
        try:
            logger.info(f"Скачиваем последовательность с Accession Number: {accession_number}")
            handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
            sequence_data = handle.read()
            handle.close()

            if not sequence_data:
                raise ValueError("Пустая последовательность")

            with open(output_file, "w") as f:
                f.write(sequence_data)

            logger.info(f"Последовательность {accession_number} сохранена в {output_file}")
            return True

        except Exception as e:
            logger.error(f"Ошибка загрузки {accession_number}: {e}")
            time.sleep(2)  # Пауза между попытками

    logger.error(f"Не удалось скачать последовательность {accession_number}")
    return False


# Функция для поиска и скачивания данных из NCBI для указанных генов
def download_sequences_for_genes(genes, organism, output_file, database="nucleotide", max_results=5, max_sequence_length=2000):
    sequences = []
    for gene in genes:
        search_query = f"{organism}[Organism] AND {gene}[Gene]"
        try:
            logger.info(f"Ищем последовательности для запроса: {search_query}")

            # Поиск записей в базе данных
            handle = Entrez.esearch(db=database, term=search_query, retmax=1000)
            record = Entrez.read(handle)
            handle.close()

            # Получение ID найденных записей
            ids = record["IdList"]
            if not ids:
                logger.warning(f"Нет найденных записей для {gene}.")
                continue

            # Получаем информацию о длине последовательностей
            handle = Entrez.esummary(db=database, id=ids)
            summaries = Entrez.read(handle)
            handle.close()

            # Фильтруем последовательности по длине
            filtered_ids = []
            for docsum in summaries['DocumentSummary']:
                seq_length = int(docsum.get('Length', 0))
                if seq_length <= max_sequence_length:
                    filtered_ids.append(docsum['Id'])
                if len(filtered_ids) >= max_results:
                    break

            if not filtered_ids:
                logger.warning(f"Нет последовательностей длиной <= {max_sequence_length} для {gene}.")
                continue

            logger.info(f"Найдено {len(filtered_ids)} последовательностей для {gene}. Скачиваем данные...")

            # Скачивание последовательностей в формате FASTA
            fasta_handle = Entrez.efetch(
                db=database, id=filtered_ids, rettype="fasta", retmode="text"
            )
            sequence_data = fasta_handle.read()
            fasta_handle.close()

            # Запись во временный файл
            temp_file = os.path.join(output_dir, f"temp_{gene}.fasta")
            with open(temp_file, "w") as f:
                f.write(sequence_data)

            # Чтение последовательностей и добавление в список
            seqs = list(SeqIO.parse(temp_file, "fasta"))
            sequences.extend(seqs)
            os.remove(temp_file)
            logger.info(f"Последовательности гена {gene} добавлены.")

        except Exception as e:
            logger.error(f"Ошибка при скачивании данных для {gene}: {e}")

    if sequences:
        # Сохраняем объединенные последовательности в output_file
        SeqIO.write(sequences, output_file, "fasta")
        logger.info(f"Все последовательности генов {genes} сохранены в {output_file}")
        return True
    else:
        logger.warning(f"Не удалось найти последовательности для генов {genes} в {organism}.")
        return False


# Функция для проверки, является ли последовательность ДНК
def is_dna(seq):
    dna_bases = set("ATCGNatcgn")
    return set(str(seq)).issubset(dna_bases)


# Функция для обрезки последовательности до длины, кратной трём
def trim_sequence(seq):
    length = len(seq)
    if length % 3 != 0:
        seq = seq[: -(length % 3)]  # Обрезаем лишние нуклеотиды
    return seq


# Функция для обработки последовательностей
def process_sequences(seq_records):
    protein_sequences = []
    for seq in seq_records:
        seq_str = str(seq.seq)
        if is_dna(seq_str):
            try:
                trimmed_seq = trim_sequence(seq.seq)
                protein_seq = trimmed_seq.translate(to_stop=True)
                protein_sequences.append(
                    SeqRecord(
                        protein_seq,
                        id=seq.id,
                        description="Translated protein sequence",
                    )
                )
            except Exception as e:
                logger.error(f"Ошибка при трансляции последовательности {seq.id}: {e}")
        else:
            protein_sequences.append(seq)
    return protein_sequences


# Функция для фильтрации последовательностей по генам
def filter_sequences_by_genes(seq_records, genes):
    filtered_records = []
    for record in seq_records:
        description = record.description.lower()
        if any(gene.lower() in description for gene in genes):
            filtered_records.append(record)
        else:
            logger.warning(f"Последовательность {record.id} не относится к генам {genes} и будет исключена.")
    return filtered_records


# Функция для получения консенсусной последовательности с использованием Bio.motifs
def get_consensus_sequence(alignment):
    # Создание объекта Motif из выравнивания
    motif = motifs.create([record.seq for record in alignment])

    # Получение консенсусной последовательности
    consensus = motif.consensus
    degenerate_consensus = motif.degenerate_consensus

    logger.info(f"Консенсусная последовательность: {consensus}")
    logger.info(f"Дегенеративная консенсусная последовательность: {degenerate_consensus}")

    return consensus, degenerate_consensus


# Функция для идентификации расширенных мутаций
def identify_advanced_mutations(sequences, reference_sequence):
    """Расширенный анализ мутаций"""
    mutations = {}
    mutation_stats = {
        'total_mutations': 0,
        'mutation_types': {
            'transitions': 0,
            'transversions': 0
        }
    }

    ref_seq = str(reference_sequence.seq)

    for record in sequences:
        if record.id == reference_sequence.id:
            continue

        seq = str(record.seq)
        muts = []

        for i, (ref_base, base) in enumerate(zip(ref_seq, seq)):
            if ref_base != base:
                muts.append((i + 1, ref_base, base))

                # Статистика типов мутаций
                mutation_stats['total_mutations'] += 1
                if (ref_base in 'AG' and base in 'AG') or (ref_base in 'CT' and base in 'CT'):
                    mutation_stats['mutation_types']['transitions'] += 1
                else:
                    mutation_stats['mutation_types']['transversions'] += 1

        mutations[record.id] = muts

    return mutations, mutation_stats


# Функция для экспорта мутаций в Excel
def export_mutations_to_excel(mutations):
    all_mutations = []
    for seq_id, seq_mutations in mutations.items():
        for pos, ref_base, base in seq_mutations:
            all_mutations.append({
                'Sequence': seq_id,
                'Position': pos,
                'Reference Base': ref_base,
                'Mutated Base': base
            })

    df = pd.DataFrame(all_mutations)
    excel_path = mutation_analysis_file
    df.to_excel(excel_path, index=False)
    logger.info(f"Мутации экспортированы в {excel_path}")


# Функция для идентификации замен аминокислот
def identify_amino_acid_changes(alignment, reference_id):
    reference_seq = None
    for record in alignment:
        if record.id == reference_id:
            reference_seq = str(record.seq)
            break
    if reference_seq is None:
        logger.error(f"Референсная белковая последовательность с ID {reference_id} не найдена.")
        return {}

    amino_acid_changes = {}
    for record in alignment:
        if record.id == reference_id:
            continue
        seq = str(record.seq)
        changes = []
        for i, (ref_aa, aa) in enumerate(zip(reference_seq, seq)):
            if ref_aa != aa and ref_aa != '-' and aa != '-':
                changes.append((i + 1, ref_aa, aa))
        amino_acid_changes[record.id] = changes
    return amino_acid_changes


# Функция для подсчета различий и схожести между группами
def compute_differences_and_similarity(ancient_seqs, vulgare_seqs, spontaneum_seqs):
    differences = {}
    similarity = {}

    group_pairs = [
        ("ancient_vs_vulgare", ancient_seqs, vulgare_seqs),
        ("ancient_vs_spontaneum", ancient_seqs, spontaneum_seqs),
        ("vulgare_vs_spontaneum", vulgare_seqs, spontaneum_seqs),
    ]

    for key, group1, group2 in group_pairs:
        diff_count = 0
        total_similarity = 0
        comparisons = 0

        for seq1 in group1:
            for seq2 in group2:
                seq1_str = str(seq1.seq)
                seq2_str = str(seq2.seq)

                # Игнорируем выравнивания, где одна из последовательностей пустая
                if not seq1_str or not seq2_str:
                    continue

                matches = sum(a == b for a, b in zip(seq1_str, seq2_str) if a != '-' and b != '-')
                alignment_length = sum(a != '-' and b != '-' for a, b in zip(seq1_str, seq2_str))

                if alignment_length > 0:
                    differences_in_pair = alignment_length - matches
                    diff_count += differences_in_pair
                    total_similarity += (matches / alignment_length) * 100
                    comparisons += 1

        differences[key] = diff_count
        similarity[key] = total_similarity / comparisons if comparisons > 0 else 0

    return differences, similarity


# Функция для реконструкции морфологии колоса и таксономической принадлежности
def reconstruct_morphology(mutations):
    with open(morphology_report_file, "w", encoding='utf-8') as f:
        f.write("Отчет о реконструкции морфологии колоса древнего ячменя\n\n")

        # Анализ мутаций в генах
        genes = ["Btr1", "Btr2", "Vrs1"]
        for gene in genes:
            f.write(f"Анализ мутаций в гене {gene}:\n")
            ancient_seq_id = f"Ancient_{gene}"
            ref_seq_id = f"{gene}_ref"
            if ref_seq_id in mutations and ancient_seq_id in mutations:
                gene_mutations = mutations[ancient_seq_id]
                if gene_mutations:
                    f.write(f"- Обнаружены мутации в гене {gene}: {len(gene_mutations)} мутаций\n")
                    for pos, ref_base, mut_base in gene_mutations:
                        f.write(f"  Позиция {pos}: {ref_base} -> {mut_base}\n")
                    # Связать мутации с фенотипическими эффектами
                    if gene in ["Btr1", "Btr2"]:
                        f.write("  Это может указывать на изменения в ломкости колоса.\n")
                    elif gene == "Vrs1":
                        f.write("  Это может влиять на форму колоса (двурядный или шестирядный).\n")
                else:
                    f.write(f"- Мутации в гене {gene} не обнаружены.\n")
            else:
                f.write(f"- Последовательности гена {gene} для анализа не найдены.\n")
            f.write("\n")

        # Общий вывод
        f.write("На основе выявленных мутаций можно предположить, что древний ячмень обладал следующими характеристиками:\n")
        f.write("- Возможные изменения в ломкости колоса, что может свидетельствовать о ранних этапах одомашнивания.\n")
        f.write("- Изменения в форме колоса, что может указывать на переход от дикого к культурному типу.\n")

        f.write("\nСравнивая с известными данными о Hordeum vulgare и Hordeum spontaneum, можно сделать вывод о таксономической принадлежности древнего ячменя.\n")
    logger.info(f"Отчет о морфологии сохранен в {morphology_report_file}")


# Функция для проведения филогенетического анализа
def perform_phylogenetic_analysis(alignment):
    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    # Сохранение дерева в изображение
    fig = plt.figure(figsize=(10, 10))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=axes)
    plt.savefig(phylogenetic_tree_file)
    logger.info(f"Филогенетическое дерево сохранено в {phylogenetic_tree_file}")

    # Анализ дерева
    ancient_ids = [record.id for record in alignment if record.id.startswith("Ancient_")]
    closest_relatives = {}
    for ancient_id in ancient_ids:
        clade = tree.find_any(name=ancient_id)
        if clade:
            # Ищем ближайших соседей среди листьев (терминальных узлов)
            distances = []
            for leaf in tree.get_terminals():
                if leaf.name != ancient_id:
                    distance = tree.distance(clade, leaf)
                    distances.append((distance, leaf.name))
            if distances:
                distances.sort()
                # Берём три ближайших родственника
                closest_relatives[ancient_id] = [name for dist, name in distances[:3]]
            else:
                closest_relatives[ancient_id] = []
        else:
            closest_relatives[ancient_id] = []

    # Получение дополнительной информации о ближайших родственниках
    relatives_info = {}
    for ancient_id, relatives in closest_relatives.items():
        info_list = []
        for relative_id in relatives:
            # Убираем префиксы, если они есть (например, "Vulgare_", "Spontaneum_")
            original_id = relative_id.split("_", 1)[-1]
            # Пытаемся получить информацию из NCBI
            try:
                # Используем esummary для получения информации
                handle = Entrez.esummary(db="nucleotide", id=original_id)
                summary = Entrez.read(handle)
                handle.close()
                if summary:
                    docsum = summary[0]
                    organism = docsum.get('Organism', 'Неизвестный организм')
                    title = docsum.get('Title', '')
                    # Таксономическая информация
                    taxid = docsum.get('TaxId', '')
                    lineage = 'Информация недоступна'
                    country = 'Информация недоступна'

                    if taxid:
                        # Получаем таксономическую информацию
                        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                        tax_records = Entrez.read(handle)
                        handle.close()
                        if tax_records:
                            lineage = tax_records[0].get('Lineage', 'Информация недоступна')
                            organism = tax_records[0].get('ScientificName', organism)

                    info_list.append({
                        "id": relative_id,
                        "organism": organism,
                        "title": title,
                        "lineage": lineage,
                        "country": country
                    })
                else:
                    info_list.append({
                        "id": relative_id,
                        "organism": "Информация недоступна",
                        "title": "Информация недоступна",
                        "lineage": "Информация недоступна",
                        "country": "Информация недоступна"
                    })
            except Exception as e:
                logger.warning(f"Не удалось получить информацию для {relative_id}: {e}")
                info_list.append({
                    "id": relative_id,
                    "organism": "Информация недоступна",
                    "title": "Информация недоступна",
                    "lineage": "Информация недоступна",
                    "country": "Информация недоступна"
                })
        relatives_info[ancient_id] = info_list

    # Записываем отчет
    with open(region_report_file, "w", encoding='utf-8') as f:
        f.write("Отчет о возможном регионе возделывания древнего ячменя\n\n")
        if relatives_info:
            for ancient_id, info_list in relatives_info.items():
                f.write(f"Для последовательности {ancient_id} ближайшими родственниками являются:\n")
                if info_list:
                    for info in info_list:
                        f.write(f"- Идентификатор: {info['id']}\n")
                        f.write(f"  Организм: {info['organism']}\n")
                        f.write(f"  Название: {info['title']}\n")
                        f.write(f"  Таксономия: {info['lineage']}\n")
                        f.write(f"  Регион: {info['country']}\n\n")
                    f.write("Это может указывать на генетическую близость и возможный регион возделывания.\n\n")
                else:
                    f.write("Ближайшие родственники не определены.\n\n")
        else:
            f.write("Не удалось определить ближайших родственников древнего ячменя на основе филогенетического анализа.\n")
    logger.info(f"Отчет о регионе возделывания сохранен в {region_report_file}")


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
    ref_accession_numbers = {
        "Btr1": "AK373702",
        "Btr2": "AK368246",
        "Vrs1": "AY093351"
    }
    reference_sequences = []
    for gene, acc_num in ref_accession_numbers.items():
        output_file = os.path.join(data_dir, f"{gene}_reference.fasta")
        if not os.path.exists(output_file):
            success = download_sequence_by_accession(acc_num, output_file)
            if not success:
                logger.error(f"Не удалось скачать референсную последовательность {gene}")
                continue

        record = SeqIO.read(output_file, "fasta")
        record.id = gene + "_ref"
        reference_sequences.append(record)

    if not reference_sequences:
        logger.error("Нет доступных референсных последовательностей. Прерывание анализа.")
        return

    # Ограничиваем количество и размер скачиваемых последовательностей
    max_results = 5
    max_sequence_length = 2000

    if not os.path.exists(vulgare_file):
        download_sequences_for_genes(genes, "Hordeum vulgare", vulgare_file, max_results=max_results, max_sequence_length=max_sequence_length)

    if not os.path.exists(spontaneum_file):
        download_sequences_for_genes(genes, "Hordeum vulgare subsp. spontaneum", spontaneum_file, max_results=max_results, max_sequence_length=max_sequence_length)

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
    spontaneum_nuc_sequences = filter_sequences_by_genes(spontaneum_nuc_sequences, genes)
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
    logger.info(f"Объединенные нуклеотидные последовательности сохранены в {combined_nuc_file}")

    # Выравнивание нуклеотидных последовательностей (Критерий 2)
    run_clustal_omega(combined_nuc_file, aligned_nuc_file, seqtype="DNA")
    alignment_nuc = AlignIO.read(aligned_nuc_file, "fasta")
    logger.info(f"Количество выравненных нуклеотидных последовательностей: {len(alignment_nuc)}")

    # Идентификация мутаций (Критерий 1)
    mutations, mutation_stats = identify_advanced_mutations(all_nuc_sequences, reference_sequences[0])
    export_mutations_to_excel(mutations)

    # Чтение и обработка белковых последовательностей
    # Трансляция нуклеотидных последовательностей в белковые
    all_protein_sequences = process_sequences(all_nuc_sequences)
    SeqIO.write(all_protein_sequences, combined_prot_file, "fasta")
    logger.info(f"Объединенные белковые последовательности сохранены в {combined_prot_file}")

    # Выравнивание белковых последовательностей (Критерий 3)
    run_clustal_omega(combined_prot_file, aligned_prot_file, seqtype="Protein")
    alignment_prot = AlignIO.read(aligned_prot_file, "fasta")
    logger.info(f"Количество выравненных белковых последовательностей: {len(alignment_prot)}")

    # Анализ замен аминокислот
    amino_acid_changes = identify_amino_acid_changes(alignment_prot, reference_sequences[0].id)
    with open(os.path.join(output_dir, "amino_acid_changes.txt"), "w", encoding='utf-8') as f:
        for seq_id, changes in amino_acid_changes.items():
            f.write(f"Замены аминокислот в последовательности {seq_id} относительно {reference_sequences[0].id}:\n")
            for pos, ref_aa, aa in changes:
                f.write(f"Позиция {pos}: {ref_aa} -> {aa}\n")
            f.write("\n")
    logger.info("Замены аминокислот сохранены в файл: amino_acid_changes.txt")

    # Получение консенсусной последовательности
    consensus, degenerate_consensus = get_consensus_sequence(alignment_prot)
    with open(consensus_file, "w", encoding='utf-8') as f:
        f.write("Консенсусная последовательность:\n")
        f.write(str(consensus) + "\n\n")
        f.write("Дегенеративная консенсусная последовательность:\n")
        f.write(str(degenerate_consensus) + "\n")
    logger.info("Консенсусные последовательности сохранены в файл: consensus_sequences.txt")

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
    ancient_seqs = [seq for seq in all_protein_sequences if seq.id.startswith("Ancient_")]
    vulgare_seqs = [seq for seq in all_protein_sequences if seq.id.startswith("Vulgare_")]
    spontaneum_seqs = [seq for seq in all_protein_sequences if seq.id.startswith("Spontaneum_")]

    # Подсчет различий и схожести
    differences, similarity = compute_differences_and_similarity(ancient_seqs, vulgare_seqs, spontaneum_seqs)

    # Сохранение результатов
    with open(differences_file, "w", encoding='utf-8') as f:
        f.write("Подсчет различий между группами:\n")
        for key, value in differences.items():
            f.write(f"{key}: {value} замен\n")
    logger.info("Таблица с подсчетом различий сохранена.")

    with open(similarity_file, "w", encoding='utf-8') as f:
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
