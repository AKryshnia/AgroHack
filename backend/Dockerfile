FROM python:3.12-slim

# Установим рабочую директорию
WORKDIR /app

# Копируем файлы проекта
COPY . /app

# Устанавливаем зависимости
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Устанавливаем Clustal Omega
RUN apt-get update && apt-get install -y clustalo && rm -rf /var/lib/apt/lists/*

# Экспонируем порт
EXPOSE 8000

# Команда запуска
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
