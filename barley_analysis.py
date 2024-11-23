from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio import AlignIO
import matplotlib.pyplot as plt
from collections import Counter
import subprocess
import os
import shutil


# Функция для проверки установки Clustal Omega
def check_clustal_omega_installed():
    if not shutil.which("clustalo"):
        raise EnvironmentError("Clustal Omega не установлен. Убедитесь, что он доступен в PATH.")


# Функция для выполнения Clustal Omega
def run_clustal_omega(input_file, output_file):
    try:
        subprocess.run(
            ["clustalo", "-i", input_file, "-o", output_file, "--auto", "--seqtype=Protein", "--force"],
            check=True,
        )
        print(f"Выравнивание завершено. Результаты сохранены в {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Ошибка при выполнении Clustal Omega: {e}")


# Функция для подсчета частоты аминокислот
def calculate_amino_acid_frequency(alignment):
    amino_acid_counts = Counter()
    for record in alignment:
        amino_acid_counts.update(str(record.seq).replace("-", ""))  # Убираем gap-символы
    total_count = sum(amino_acid_counts.values())
    frequencies = {aa: count / total_count for aa, count in amino_acid_counts.items()}
    return frequencies


# Главная функция
def main():
    # Пути к файлам
    ancient_barley_nuc_file = "data/ancient_barley.fasta"
    vulgare_file = "data/protein.faa"
    spontaneum_file = "data/spontaneum_genes.fasta"
    combined_file = "combined_sequences.fasta"
    aligned_file = "aligned_sequences.fasta"

    # Проверяем установку Clustal Omega
    check_clustal_omega_installed()

    # Чтение древних последовательностей и трансляция в белковые
    ancient_nuc_sequences = list(SeqIO.parse(ancient_barley_nuc_file, "fasta"))
    ancient_protein_sequences = [
        SeqRecord(seq.seq.translate(to_stop=True), id=seq.id, description="Translated protein sequence")
        for seq in ancient_nuc_sequences
    ]
    print(f"Найдено последовательностей в древнем ячмене: {len(ancient_protein_sequences)}")
    for i, protein in enumerate(ancient_protein_sequences):
        print(f"Длина белковой последовательности {i + 1}: {len(protein.seq)}")

    # Чтение современных последовательностей
    vulgare_sequences = list(SeqIO.parse(vulgare_file, "fasta"))
    spontaneum_sequences = list(SeqIO.parse(spontaneum_file, "fasta"))

    # Объединение всех последовательностей
    all_sequences = ancient_protein_sequences + vulgare_sequences + spontaneum_sequences
    SeqIO.write(all_sequences, combined_file, "fasta")
    print(f"Объединенные последовательности сохранены в {combined_file}")

    # Запуск Clustal Omega для выравнивания
    run_clustal_omega(combined_file, aligned_file)

    # Чтение выравнивания
    alignment = AlignIO.read(aligned_file, "fasta")
    print(f"Количество выравненных последовательностей: {len(alignment)}")
    print(f"Длина выравнивания: {alignment.get_alignment_length()}")

    # Получение консенсусной последовательности
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    print(f"Консенсусная последовательность: {consensus}")

    # Частота аминокислот
    frequencies = calculate_amino_acid_frequency(alignment)

    # Построение гистограммы
    amino_acids = list(frequencies.keys())
    freq_values = list(frequencies.values())
    plt.figure(figsize=(10, 6))
    plt.bar(amino_acids, freq_values)
    plt.xlabel("Аминокислоты")
    plt.ylabel("Частота")
    plt.title("Частота аминокислот в выравненных последовательностях")
    plt.savefig("amino_acid_frequencies.png")
    print("График сохранён в файл: amino_acid_frequencies.png")


if __name__ == "__main__":
    main()
