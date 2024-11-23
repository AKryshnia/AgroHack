const API_BASE_URL = import.meta.env.VITE_API_BASE_URL;

export const uploadFile = async (file) => {
  const formData = new FormData();
  formData.append('file', file);

  const response = await fetch(`${API_BASE_URL}api/upload-csv/`, {
    method: 'POST',
    body: formData,
  });

  if (!response.ok) {
    const errorText = await response.text();
    throw new Error(errorText || 'Ошибка при загрузке файла');
  }

  const data = await response.json();
  return data;
};
