# barley_protein_analysis.py

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
from collections import Counter
import os
import shutil

# Пути к файлам
ancient_barley_nuc_file = 'ancient_barley.fasta'  # Нуклеотидная последовательность древнего ячменя
hordeum_vulgare_protein_file = 'data/protein.faa'  # Белковые последовательности Hordeum vulgare
hordeum_spontaneum_protein_file = 'spontaneum_protein.faa'  # Если есть данные

# Целевые гены для фильтрации
target_genes = ['Btr1', 'Btr2', 'Vrs1']


def read_nucleotide_sequence(file_path):
    sequence = SeqIO.read(file_path, 'fasta')
    return sequence


def translate_sequence(nuc_record):
    # Используем стандартный генетический код (table=1)
    protein_seq = nuc_record.seq.translate(to_stop=True)
    protein_record = SeqRecord(protein_seq, id=nuc_record.id, description='Translated protein sequence')
    return protein_record


def read_protein_sequences(file_path):
    sequences = list(SeqIO.parse(file_path, 'fasta'))
    return sequences


def filter_sequences_by_genes(sequences, target_genes):
    filtered_sequences = []
    for record in sequences:
        for gene in target_genes:
            if gene.lower() in record.description.lower():
                filtered_sequences.append(record)
                break
    return filtered_sequences


def calculate_amino_acid_frequency(alignment):
    amino_acid_counts = Counter()
    for record in alignment:
        amino_acid_counts.update(str(record.seq).replace('-', ''))  # Убираем gap-символы
    total_count = sum(amino_acid_counts.values())
    frequencies = {aa: count / total_count for aa, count in amino_acid_counts.items()}
    return frequencies


def main():
    # Чтение и трансляция древней последовательности
    ancient_nuc_sequence = read_nucleotide_sequence(ancient_barley_nuc_file)
    ancient_protein_sequence = translate_sequence(ancient_nuc_sequence)
    print(f'Длина аминокислотной последовательности древнего ячменя: {len(ancient_protein_sequence.seq)}')
    
    # Чтение и фильтрация белковых последовательностей Hordeum vulgare
    vulgare_protein_sequences = read_protein_sequences(hordeum_vulgare_protein_file)
    print(f'Общее количество белковых последовательностей Hordeum vulgare: {len(vulgare_protein_sequences)}')
    filtered_vulgare_proteins = filter_sequences_by_genes(vulgare_protein_sequences, target_genes)
    print(f'Количество белковых последовательностей Hordeum vulgare после фильтрации: {len(filtered_vulgare_proteins)}')
    
    # Чтение и фильтрация белковых последовательностей Hordeum spontaneum (если доступно)
    if os.path.exists(hordeum_spontaneum_protein_file):
        spontaneum_protein_sequences = read_protein_sequences(hordeum_spontaneum_protein_file)
        filtered_spontaneum_proteins = filter_sequences_by_genes(spontaneum_protein_sequences, target_genes)
        print(
            f'Количество белковых последовательностей Hordeum spontaneum после фильтрации: {len(filtered_spontaneum_proteins)}')
    else:
        filtered_spontaneum_proteins = []
        print('Файл с белковыми последовательностями Hordeum spontaneum не найден. Пропускаем этот шаг.')
    
    # Собираем все последовательности для выравнивания
    sequences_to_align = [ancient_protein_sequence] + filtered_vulgare_proteins + filtered_spontaneum_proteins
    
    if len(sequences_to_align) < 2:
        print('Недостаточно последовательностей для выравнивания. Завершение программы.')
        return
    
    # Сохраняем последовательности в файл
    combined_protein_sequences_file = 'combined_protein_sequences.fasta'
    SeqIO.write(sequences_to_align, combined_protein_sequences_file, 'fasta')
    print(f'Объединенные аминокислотные последовательности сохранены в {combined_protein_sequences_file}')
    
    # Проверяем наличие Clustal Omega
    clustalo_path = 'clustalo'
    if not shutil.which(clustalo_path):
        print('Clustal Omega не найден. Убедитесь, что он установлен и доступен в PATH.')
        return
    
    # Выполняем выравнивание
    aligned_protein_sequences_file = 'aligned_protein_sequences.fasta'
    clustalomega_cline = ClustalOmegaCommandline(cmd=clustalo_path, infile=combined_protein_sequences_file,
                                                 outfile=aligned_protein_sequences_file, verbose=True, auto=True,
                                                 seqtype="Protein")
    stdout, stderr = clustalomega_cline()
    print(f'Выравнивание завершено. Результаты сохранены в {aligned_protein_sequences_file}')
    
    # Чтение выравнивания
    alignment = AlignIO.read(aligned_protein_sequences_file, 'fasta')
    print(f'Количество выравненных последовательностей: {len(alignment)}')
    print(f'Длина выравнивания: {alignment.get_alignment_length()}')
    
    # Получаем консенсусную последовательность
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    print('Консенсусная последовательность:')
    print(consensus)
    
    # Вычисляем частоты аминокислот
    frequencies = calculate_amino_acid_frequency(alignment)
    
    # Построение гистограммы
    amino_acids = list(frequencies.keys())
    freq_values = list(frequencies.values())
    
    plt.figure(figsize=(10, 6))
    plt.bar(amino_acids, freq_values)
    plt.xlabel('Аминокислоты')
    plt.ylabel('Частота')
    plt.title('Частота аминокислот в выравненных последовательностях')
    plt.show()


if __name__ == '__main__':
    main()
