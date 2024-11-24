# Makefile для Ruff

# Настройки
PYTHON_FILES = .  # Укажите директорию или файлы для проверки
RUFF_ARGS = --fix  # Аргументы для Ruff (можно убрать --fix для проверки без исправлений)

# Цели
.PHONY: help install lint fix clean

help:
	@echo "Доступные команды:"
	@echo "  make install  - Установить зависимости для Ruff"
	@echo "  make lint     - Проверить код с помощью Ruff"
	@echo "  make fix      - Исправить проблемы линтинга автоматически"
	@echo "  make clean    - Очистить временные файлы (например, .pyc)"

install:
	@echo "Устанавливаем Ruff..."
	pip install ruff

lint:
	@echo "Запускаем линтинг..."
	ruff $(PYTHON_FILES)

fix:
	@echo "Исправляем код..."
	ruff $(PYTHON_FILES) $(RUFF_ARGS)

clean:
	@echo "Очищаем временные файлы..."
	find . -name "*.pyc" -o -name "__pycache__" | xargs rm -rf
