# Stage 1: build
FROM node:18-alpine as build

WORKDIR /app

# Копируем package.json и package-lock.json
COPY package.json package-lock.json ./

# Устанавливаем зависимости
RUN npm install

# Копируем остальные файлы
COPY . .

# Сборка приложения
RUN npm run build

# Stage 2: serve
FROM nginx:stable-alpine

# Копируем собранные файлы из Stage 1
COPY --from=build /app/dist /usr/share/nginx/html

# Копируем конфигурацию Nginx
COPY nginx.conf /etc/nginx/conf.d/default.conf

# Экспонируем порт
EXPOSE 80

# Команда запуска
CMD ["nginx", "-g", "daemon off;"]
