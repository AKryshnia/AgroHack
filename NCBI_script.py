from Bio import Entrez
import os

# Укажите свою почту для работы с NCBI
Entrez.email = "your_email@example.com"


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


# Пример использования
if __name__ == "__main__":
    # Создаем папку data, если ее нет
    if not os.path.exists("data"):
        os.makedirs("data")

    # Запрос: Hordeum vulgare subsp. spontaneum и гены Btr1, Btr2, Vrs1
    query_spontaneum = (
        "Hordeum vulgare subsp. spontaneum[Organism] AND (Btr1 OR Btr2 OR Vrs1)"
    )
    output_spontaneum = "data/spontaneum_genes.fasta"
    download_sequences(
        search_query=query_spontaneum, output_file=output_spontaneum, max_results=20
    )

    # Запрос: Hordeum vulgare и гены Btr1, Btr2, Vrs1
    query_vulgare = "Hordeum vulgare[Organism] AND (Btr1 OR Btr2 OR Vrs1)"
    output_vulgare = "data/hordeum_vulgare.fasta"
    download_sequences(
        search_query=query_vulgare, output_file=output_vulgare, max_results=20
    )
