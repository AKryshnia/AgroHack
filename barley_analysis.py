from Bio import SeqIO, Entrez, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio import motifs
import matplotlib.pyplot as plt
from collections import Counter
import subprocess
import os
import shutil

# Укажите свою почту для работы с NCBI
Entrez.email = "your_email@example.com"


# Функция для проверки установки Clustal Omega
def check_clustal_omega_installed():
    if not shutil.which("clustalo"):
        raise EnvironmentError(
            "Clustal Omega не установлен. Убедитесь, что он доступен в PATH."
        )


# Функция для выполнения Clustal Omega
def run_clustal_omega(input_file, output_file):
    try:
        subprocess.run(
            [
                "clustalo",
                "-i",
                input_file,
                "-o",
                output_file,
                "--auto",
                "--seqtype=Protein",
                "--force",
            ],
            check=True,
        )
        print(f"Выравнивание завершено. Результаты сохранены в {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Ошибка при выполнении Clustal Omega: {e}")


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


# Функция для поиска и скачивания данных из NCBI
def download_sequences(
    search_query, output_file, database="nucleotide", max_results=10
):
    try:
        print(f"Ищем последовательности для запроса: {search_query}")

        # Поиск записей в базе данных
        handle = Entrez.esearch(db=database, term=search_query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        # Получение ID найденных записей
        ids = record["IdList"]
        if not ids:
            print("Нет найденных записей. Проверьте запрос.")
            return

        print(f"Найдено {len(ids)} записей. Скачиваем данные...")

        # Скачивание последовательностей в формате FASTA
        fasta_handle = Entrez.efetch(
            db=database, id=ids, rettype="fasta", retmode="text"
        )
        with open(output_file, "w") as f:
            f.write(fasta_handle.read())
        fasta_handle.close()

        print(f"Данные успешно сохранены в {output_file}")
    except Exception as e:
        print(f"Ошибка при скачивании данных: {e}")


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
                print(f"Ошибка при трансляции последовательности {seq.id}: {e}")
        else:
            protein_sequences.append(seq)
    return protein_sequences


# Функция для получения консенсусной последовательности с использованием Bio.motifs
def get_consensus_sequence(alignment):
    # Создание объекта Motif из выравнивания
    motif = motifs.create([record.seq for record in alignment])

    # Получение консенсусной последовательности
    consensus = motif.consensus
    degenerate_consensus = motif.degenerate_consensus  # Убрали вызов как метода

    print(f"Консенсусная последовательность: {consensus}")
    print(f"Дегенеративная консенсусная последовательность: {degenerate_consensus}")

    return consensus, degenerate_consensus


# Главная функция
def main():
    # Проверяем наличие директории data
    if not os.path.exists("data"):
        os.makedirs("data")

    # Пути к файлам
    ancient_barley_nuc_file = "data/ancient_barley.fasta"
    vulgare_file = "data/hordeum_vulgare.fasta"
    spontaneum_file = "data/spontaneum_genes.fasta"
    combined_file = "combined_sequences.fasta"
    aligned_file = "aligned_sequences.fasta"

    # Проверяем установку Clustal Omega
    check_clustal_omega_installed()

    # Скачиваем данные, если файлы не существуют
    if not os.path.exists(vulgare_file):
        # Запросы для современного ячменя
        query_vulgare = "Hordeum vulgare[Organism] AND (Btr1 OR Btr2 OR Vrs1)"
        download_sequences(
            search_query=query_vulgare, output_file=vulgare_file, max_results=20
        )

    if not os.path.exists(spontaneum_file):
        # Запросы для дикого ячменя
        query_spontaneum = (
            "Hordeum vulgare subsp. spontaneum[Organism] AND (Btr1 OR Btr2 OR Vrs1)"
        )
        download_sequences(
            search_query=query_spontaneum, output_file=spontaneum_file, max_results=20
        )

    # Чтение древних последовательностей и трансляция в белковые
    if not os.path.exists(ancient_barley_nuc_file):
        print(
            f"Файл {ancient_barley_nuc_file} не найден. Пожалуйста, поместите его в директорию data."
        )
        return

    ancient_nuc_sequences = list(SeqIO.parse(ancient_barley_nuc_file, "fasta"))
    ancient_protein_sequences = process_sequences(ancient_nuc_sequences)
    print(
        f"Найдено последовательностей в древнем ячмене: {len(ancient_protein_sequences)}"
    )
    for i, protein in enumerate(ancient_protein_sequences):
        print(f"Длина белковой последовательности {i + 1}: {len(protein.seq)}")

    # Чтение современных последовательностей и трансляция при необходимости
    vulgare_nuc_sequences = list(SeqIO.parse(vulgare_file, "fasta"))
    vulgare_protein_sequences = process_sequences(vulgare_nuc_sequences)

    spontaneum_nuc_sequences = list(SeqIO.parse(spontaneum_file, "fasta"))
    spontaneum_protein_sequences = process_sequences(spontaneum_nuc_sequences)

    # Объединение всех последовательностей
    all_sequences = (
        ancient_protein_sequences
        + vulgare_protein_sequences
        + spontaneum_protein_sequences
    )
    SeqIO.write(all_sequences, combined_file, "fasta")
    print(f"Объединенные последовательности сохранены в {combined_file}")

    # Запуск Clustal Omega для выравнивания
    run_clustal_omega(combined_file, aligned_file)

    # Чтение выравнивания
    alignment = AlignIO.read(aligned_file, "fasta")
    print(f"Количество выравненных последовательностей: {len(alignment)}")
    print(f"Длина выравнивания: {alignment.get_alignment_length()}")

    # Получение консенсусной последовательности
    consensus, degenerate_consensus = get_consensus_sequence(alignment)

    # Частота аминокислот
    frequencies = calculate_amino_acid_frequency(alignment)

    # Построение гистограммы
    amino_acids = sorted(frequencies.keys())
    freq_values = [frequencies[aa] for aa in amino_acids]
    plt.figure(figsize=(12, 6))
    plt.bar(amino_acids, freq_values, color="skyblue")
    plt.xlabel("Аминокислоты")
    plt.ylabel("Частота")
    plt.title("Частота аминокислот в выравненных последовательностях")
    plt.tight_layout()
    plt.savefig("amino_acid_frequencies.png")
    print("График сохранён в файл: amino_acid_frequencies.png")

    # Сохранение консенсусных последовательностей в файл
    with open("consensus_sequences.txt", "w") as f:
        f.write("Консенсусная последовательность:\n")
        f.write(str(consensus) + "\n\n")
        f.write("Дегенеративная консенсусная последовательность:\n")
        f.write(str(degenerate_consensus) + "\n")
    print("Консенсусные последовательности сохранены в файл: consensus_sequences.txt")

    # Подсчет различий между группами
    # Предполагается, что последовательности разделены по подспособностям в ID
    differences = {
        "ancient_vs_vulgare": 0,
        "ancient_vs_spontaneum": 0,
        "vulgare_vs_spontaneum": 0,
    }

    # Разделение последовательностей по группам
    ancient_seqs = ancient_protein_sequences
    vulgare_seqs = vulgare_protein_sequences
    spontaneum_seqs = spontaneum_protein_sequences

    # Функция для подсчета различий между двумя группами
    def count_differences(group1, group2):
        count = 0
        for seq1 in group1:
            for seq2 in group2:
                for a, b in zip(seq1.seq, seq2.seq):
                    if a != b:
                        count += 1
        return count

    differences["ancient_vs_vulgare"] = count_differences(ancient_seqs, vulgare_seqs)
    differences["ancient_vs_spontaneum"] = count_differences(
        ancient_seqs, spontaneum_seqs
    )
    differences["vulgare_vs_spontaneum"] = count_differences(
        vulgare_seqs, spontaneum_seqs
    )

    # Сохранение таблицы с подсчетом различий
    with open("differences_table.txt", "w") as f:
        f.write("Подсчет различий между группами:\n")
        for key, value in differences.items():
            f.write(f"{key}: {value} замен\n")
    print("Таблица с подсчетом различий сохранена в файл: differences_table.txt")


if __name__ == "__main__":
    main()
