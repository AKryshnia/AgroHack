# main.py

from fastapi import FastAPI, File, UploadFile, HTTPException
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
import subprocess
from pathlib import Path
import shutil

app = FastAPI()

# Разрешенные источники для CORS
origins = [
    "http://localhost",
    "http://localhost:3000",
    # Добавьте другие, если необходимо
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,  # Разрешенные источники
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Путь к директории данных и результатам
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "output"

# Убедимся, что директории существуют
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Функция для запуска barley_analysis.py
def run_analysis():
    try:
        result = subprocess.run(
            ["python", "barley_analysis.py"],
            cwd=BASE_DIR,
            capture_output=True,
            text=True,
            check=True,
        )
        print(result.stdout)
        if result.stderr:
            print("stderr:", result.stderr)
    except subprocess.CalledProcessError as e:
        print("Ошибка при выполнении barley_analysis.py:", e)
        print("stdout:", e.stdout)
        print("stderr:", e.stderr)
        raise HTTPException(status_code=500, detail="Ошибка при выполнении анализа")


@app.post("/api/run-analysis/")
async def upload_file(file: UploadFile = File(...)):
    # Проверяем тип файла (ожидаем .fasta)
    if not (
        file.content_type
        in ["application/fasta", "application/octet-stream", "text/plain"]
        or file.filename.endswith(".fasta")
    ):
        raise HTTPException(
            status_code=400,
            detail="Неверный формат файла. Пожалуйста, загрузите файл в формате .fasta",
        )

    # Сохраняем загруженный файл как ancient_barley.fasta
    try:
        ancient_barley_path = DATA_DIR / "ancient_barley.fasta"
        with ancient_barley_path.open("wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Не удалось сохранить файл: {e}")

    # Запуск анализа
    run_analysis()

    # Собираем список файлов из директории output
    generated_files = [f.name for f in OUTPUT_DIR.iterdir() if f.is_file()]

    # Вернём список файлов для скачивания
    return {"message": "Анализ выполнен успешно", "files": generated_files}


@app.get("/api/download/{filename}")
def download_file(filename: str):
    file_path = OUTPUT_DIR / filename

    if not file_path.exists():
        raise HTTPException(status_code=404, detail="Файл не найден")

    return FileResponse(
        path=file_path, filename=filename, media_type="application/octet-stream"
    )
