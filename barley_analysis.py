from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
from collections import Counter
import subprocess
import shutil  # Добавлен импорт shutil
import os

# Пути к файлам
ancient_barley_nuc_file = 'data/ancient_barley.fasta'
hordeum_vulgare_protein_file = 'data/protein.faa'
hordeum_spontaneum_protein_file = 'data/spontaneum_protein.faa'

# Целевые гены
target_genes = ['Btr1', 'Btr2', 'Vrs1']


def clean_sequence(sequence):
    """Очищает последовательность от недопустимых символов"""
    valid_chars = set("ATGCatgc")
    return "".join([char for char in sequence if char in valid_chars])


def read_nucleotide_sequences(file_path):
    """Чтение всех последовательностей из файла FASTA и очистка"""
    sequences = []
    for record in SeqIO.parse(file_path, 'fasta'):
        cleaned_seq = clean_sequence(str(record.seq))
        record.seq = cleaned_seq  # Замена на очищенную последовательность
        sequences.append(record)
    return sequences


def translate_sequence(nuc_record):
    """Трансляция нуклеотидной последовательности в белковую"""
    protein_seq = nuc_record.seq.translate(to_stop=True)
    protein_record = SeqRecord(protein_seq, id=nuc_record.id, description='Translated protein sequence')
    return protein_record


def check_clustal_omega_installed():
    """Проверяет, установлен ли Clustal Omega"""
    if not shutil.which('clustalo'):
        raise EnvironmentError("Clustal Omega не найден. Убедитесь, что он установлен и доступен в PATH.")


def run_clustal_omega(input_file, output_file):
    """Выполняет выравнивание с помощью Clustal Omega"""
    try:
        subprocess.run(
            ['clustalo', '-i', input_file, '-o', output_file, '--auto', '--seqtype=Protein'],
            check=True
        )
        print(f'Выравнивание завершено. Результаты сохранены в {output_file}')
    except subprocess.CalledProcessError as e:
        print(f'Ошибка при выполнении Clustal Omega: {e}')


def calculate_amino_acid_frequency(alignment_file):
    """Расчёт частоты аминокислот в выравненных последовательностях"""
    alignment = SeqIO.parse(alignment_file, 'fasta')
    amino_acid_counts = Counter()
    for record in alignment:
        amino_acid_counts.update(str(record.seq).replace('-', ''))
    total_count = sum(amino_acid_counts.values())
    frequencies = {aa: count / total_count for aa, count in amino_acid_counts.items()}
    return frequencies


def plot_amino_acid_frequencies(frequencies):
    """Визуализация частоты аминокислот"""
    amino_acids = list(frequencies.keys())
    freq_values = list(frequencies.values())
    plt.figure(figsize=(10, 6))
    plt.bar(amino_acids, freq_values)
    plt.xlabel('Аминокислоты')
    plt.ylabel('Частота')
    plt.title('Частота аминокислот в выравненных последовательностях')
    plt.show()


def main():
    # Проверка Clustal Omega
    check_clustal_omega_installed()

    # Чтение и трансляция древних последовательностей
    ancient_nuc_sequences = read_nucleotide_sequences(ancient_barley_nuc_file)
    print(f'Найдено последовательностей в древнем ячмене: {len(ancient_nuc_sequences)}')

    ancient_protein_sequences = [
        translate_sequence(nuc_record) for nuc_record in ancient_nuc_sequences
    ]
    for i, protein in enumerate(ancient_protein_sequences):
        print(f'Длина белковой последовательности {i + 1}: {len(protein.seq)}')

    # Сохранение белковых последовательностей в файл
    combined_file = 'combined_sequences.fasta'
    SeqIO.write(ancient_protein_sequences, combined_file, 'fasta')

    # Выполнение выравнивания
    aligned_file = 'aligned_sequences.fasta'
    run_clustal_omega(combined_file, aligned_file)

    # Анализ частоты аминокислот
    frequencies = calculate_amino_acid_frequency(aligned_file)

    # Построение гистограммы
    plot_amino_acid_frequencies(frequencies)


if __name__ == '__main__':
    main()
