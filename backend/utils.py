import os
import sys
import shutil
import logging
import subprocess
import time
from collections import Counter
from Bio import SeqIO, Entrez, Phylo
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

from folders import (
    mutation_analysis_file,
    morphology_report_file,
    phylogenetic_tree_file,
    region_report_file,
)


# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s: %(message)s",
    filename="analysis.log",
)
logger = logging.getLogger(__name__)


def create_dir_if_not_exists(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
        logger.info(f"Директория {dir} была создана.")
    else:
        logger.info(f"Директория {dir} уже существует.")


# Функция для проверки установки Clustal Omega
def check_clustal_omega_installed():
    if not shutil.which("clustalo"):
        raise EnvironmentError(
            "Clustal Omega не установлен. Убедитесь, что он доступен в PATH."
        )


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
            logger.info(
                f"Скачиваем последовательность с Accession Number: {accession_number}"
            )
            handle = Entrez.efetch(
                db="nucleotide", id=accession_number, rettype="fasta", retmode="text"
            )
            sequence_data = handle.read()
            handle.close()

            if not sequence_data:
                raise ValueError("Пустая последовательность")

            with open(output_file, "w") as f:
                f.write(sequence_data)

            logger.info(
                f"Последовательность {accession_number} сохранена в {output_file}"
            )
            return True

        except Exception as e:
            logger.error(f"Ошибка загрузки {accession_number}: {e}")
            time.sleep(2)  # Пауза между попытками

    logger.error(f"Не удалось скачать последовательность {accession_number}")
    return False


def fetch_gene_sequences(gene, organism, database, max_results, max_sequence_length):
    """Скачивает последовательности для одного гена."""
    try:
        search_query = f"{organism}[Organism] AND {gene}[Gene]"
        logger.info(f"Ищем последовательности для запроса: {search_query}")

        # Поиск записей в базе данных
        handle = Entrez.esearch(db=database, term=search_query, retmax=1000)
        record = Entrez.read(handle)
        handle.close()

        # Получение ID найденных записей
        ids = record["IdList"]
        if not ids:
            logger.warning(f"Нет найденных записей для {gene}.")
            return []

        # Получаем информацию о длине последовательностей
        handle = Entrez.esummary(db=database, id=ids)
        summaries = Entrez.read(handle)
        handle.close()

        # Фильтруем последовательности по длине
        filtered_ids = []
        for docsum in summaries["DocumentSummary"]:
            seq_length = int(docsum.get("Length", 0))
            if seq_length <= max_sequence_length:
                filtered_ids.append(docsum["Id"])
            if len(filtered_ids) >= max_results:
                break

        if not filtered_ids:
            logger.warning(
                f"Нет последовательностей длиной <= {max_sequence_length} для {gene}."
            )
            return []

        logger.info(
            f"Найдено {len(filtered_ids)} последовательностей для {gene}. Скачиваем данные..."
        )

        # Скачивание последовательностей в формате FASTA
        fasta_handle = Entrez.efetch(
            db=database, id=filtered_ids, rettype="fasta", retmode="text"
        )
        sequence_data = fasta_handle.read()
        fasta_handle.close()

        # Чтение последовательностей в список
        seqs = list(SeqIO.parse(sequence_data.splitlines(), "fasta"))
        logger.info(f"Последовательности гена {gene} добавлены.")
        return seqs

    except Exception as e:
        logger.error(f"Ошибка при скачивании данных для {gene}: {e}")
        return []


def download_sequences_for_genes(
    genes,
    organism,
    output_file,
    database="nucleotide",
    max_results=5,
    max_sequence_length=2000,
    max_workers=5,
):
    """Скачивает последовательности для списка генов с использованием многопоточности."""
    sequences = []

    # Используем ThreadPoolExecutor для параллельной загрузки данных
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_gene = {
            executor.submit(
                fetch_gene_sequences,
                gene,
                organism,
                database,
                max_results,
                max_sequence_length,
            ): gene
            for gene in genes
        }

        for future in as_completed(future_to_gene):
            gene = future_to_gene[future]
            try:
                seqs = future.result()
                sequences.extend(seqs)
            except Exception as e:
                logger.error(f"Ошибка при обработке гена {gene}: {e}")

    if sequences:
        # Сохраняем объединенные последовательности в output_file
        SeqIO.write(sequences, output_file, "fasta")
        logger.info(f"Все последовательности генов {genes} сохранены в {output_file}")
        return True
    else:
        logger.warning(
            f"Не удалось найти последовательности для генов {genes} в {organism}."
        )
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
            logger.warning(
                f"Последовательность {record.id} не относится к генам {genes} и будет исключена."
            )
    return filtered_records


# Функция для получения консенсусной последовательности с использованием Bio.motifs
def get_consensus_sequence(alignment):
    # Создание объекта Motif из выравнивания
    motif = motifs.create([record.seq for record in alignment])

    # Получение консенсусной последовательности
    consensus = motif.consensus
    degenerate_consensus = motif.degenerate_consensus

    logger.info(f"Консенсусная последовательность: {consensus}")
    logger.info(
        f"Дегенеративная консенсусная последовательность: {degenerate_consensus}"
    )

    return consensus, degenerate_consensus


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


# Функция для идентификации расширенных мутаций
def identify_advanced_mutations(sequences, reference_sequence):
    """Расширенный анализ мутаций"""
    mutations = {}
    mutation_stats = {
        "total_mutations": 0,
        "mutation_types": {"transitions": 0, "transversions": 0},
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
                mutation_stats["total_mutations"] += 1
                if (ref_base in "AG" and base in "AG") or (
                    ref_base in "CT" and base in "CT"
                ):
                    mutation_stats["mutation_types"]["transitions"] += 1
                else:
                    mutation_stats["mutation_types"]["transversions"] += 1

        mutations[record.id] = muts

    return mutations, mutation_stats


# Функция для экспорта мутаций в Excel
def export_mutations_to_excel(mutations):
    all_mutations = []
    for seq_id, seq_mutations in mutations.items():
        for pos, ref_base, base in seq_mutations:
            all_mutations.append(
                {
                    "Sequence": seq_id,
                    "Position": pos,
                    "Reference Base": ref_base,
                    "Mutated Base": base,
                }
            )

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
        logger.error(
            f"Референсная белковая последовательность с ID {reference_id} не найдена."
        )
        return {}

    amino_acid_changes = {}
    for record in alignment:
        if record.id == reference_id:
            continue
        seq = str(record.seq)
        changes = []
        for i, (ref_aa, aa) in enumerate(zip(reference_seq, seq)):
            if ref_aa != aa and ref_aa != "-" and aa != "-":
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

                matches = sum(
                    a == b for a, b in zip(seq1_str, seq2_str) if a != "-" and b != "-"
                )
                alignment_length = sum(
                    a != "-" and b != "-" for a, b in zip(seq1_str, seq2_str)
                )

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
    with open(morphology_report_file, "w", encoding="utf-8") as f:
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
                    f.write(
                        f"- Обнаружены мутации в гене {gene}: {len(gene_mutations)} мутаций\n"
                    )
                    for pos, ref_base, mut_base in gene_mutations:
                        f.write(f"  Позиция {pos}: {ref_base} -> {mut_base}\n")
                    # Связать мутации с фенотипическими эффектами
                    if gene in ["Btr1", "Btr2"]:
                        f.write(
                            "  Это может указывать на изменения в ломкости колоса.\n"
                        )
                    elif gene == "Vrs1":
                        f.write(
                            "  Это может влиять на форму колоса (двурядный или шестирядный).\n"
                        )
                else:
                    f.write(f"- Мутации в гене {gene} не обнаружены.\n")
            else:
                f.write(f"- Последовательности гена {gene} для анализа не найдены.\n")
            f.write("\n")

        # Общий вывод
        f.write(
            "На основе выявленных мутаций можно предположить, что древний ячмень обладал следующими характеристиками:\n"
        )
        f.write(
            "- Возможные изменения в ломкости колоса, что может свидетельствовать о ранних этапах одомашнивания.\n"
        )
        f.write(
            "- Изменения в форме колоса, что может указывать на переход от дикого к культурному типу.\n"
        )

        f.write(
            "\nСравнивая с известными данными о Hordeum vulgare и Hordeum spontaneum, можно сделать вывод о таксономической принадлежности древнего ячменя.\n"
        )
    logger.info(f"Отчет о морфологии сохранен в {morphology_report_file}")


# Функция для проведения филогенетического анализа
def perform_phylogenetic_analysis(alignment):
    calculator = DistanceCalculator("blosum62")
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
    ancient_ids = [
        record.id for record in alignment if record.id.startswith("Ancient_")
    ]
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
                    organism = docsum.get("Organism", "Неизвестный организм")
                    title = docsum.get("Title", "")
                    # Таксономическая информация
                    taxid = docsum.get("TaxId", "")
                    lineage = "Информация недоступна"
                    country = "Информация недоступна"

                    if taxid:
                        # Получаем таксономическую информацию
                        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                        tax_records = Entrez.read(handle)
                        handle.close()
                        if tax_records:
                            lineage = tax_records[0].get(
                                "Lineage", "Информация недоступна"
                            )
                            organism = tax_records[0].get("ScientificName", organism)

                    info_list.append(
                        {
                            "id": relative_id,
                            "organism": organism,
                            "title": title,
                            "lineage": lineage,
                            "country": country,
                        }
                    )
                else:
                    info_list.append(
                        {
                            "id": relative_id,
                            "organism": "Информация недоступна",
                            "title": "Информация недоступна",
                            "lineage": "Информация недоступна",
                            "country": "Информация недоступна",
                        }
                    )
            except Exception as e:
                logger.warning(f"Не удалось получить информацию для {relative_id}: {e}")
                info_list.append(
                    {
                        "id": relative_id,
                        "organism": "Информация недоступна",
                        "title": "Информация недоступна",
                        "lineage": "Информация недоступна",
                        "country": "Информация недоступна",
                    }
                )
        relatives_info[ancient_id] = info_list

    # Записываем отчет
    with open(region_report_file, "w", encoding="utf-8") as f:
        f.write("Отчет о возможном регионе возделывания древнего ячменя\n\n")
        if relatives_info:
            for ancient_id, info_list in relatives_info.items():
                f.write(
                    f"Для последовательности {ancient_id} ближайшими родственниками являются:\n"
                )
                if info_list:
                    for info in info_list:
                        f.write(f"- Идентификатор: {info['id']}\n")
                        f.write(f"  Организм: {info['organism']}\n")
                        f.write(f"  Название: {info['title']}\n")
                        f.write(f"  Таксономия: {info['lineage']}\n")
                        f.write(f"  Регион: {info['country']}\n\n")
                    f.write(
                        "Это может указывать на генетическую близость и возможный регион возделывания.\n\n"
                    )
                else:
                    f.write("Ближайшие родственники не определены.\n\n")
        else:
            f.write(
                "Не удалось определить ближайших родственников древнего ячменя на основе филогенетического анализа.\n"
            )
    logger.info(f"Отчет о регионе возделывания сохранен в {region_report_file}")
