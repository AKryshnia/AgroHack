import React, { useState } from 'react';
import { toast, ToastContainer } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';
import { uploadFile } from '../api';

const UploadForm = () => {
  const [file, setFile] = useState(null);
  const [downloadUrl, setDownloadUrl] = useState('');
  const [isUploaded, setIsUploaded] = useState(false);
  const [isLoading, setIsLoading] = useState(false);

  const handleFileChange = (event) => {
    const selectedFile = event.target.files[0];

    if (!selectedFile) return;

    if (selectedFile.type !== 'text/csv' && !selectedFile.name.endsWith('.csv')) {
      toast.error('Неверный формат файла. Пожалуйста, выберите файл в формате .csv.');
      return;
    }

    setFile(selectedFile);
  };

  const handleUpload = async () => {
    if (!file) {
      toast.error('Пожалуйста, выберите файл для загрузки.');
      return;
    }

    setIsLoading(true);

    try {
      const data = await uploadFile(file);

      toast.success('Файл успешно загружен и обработан.');

      setIsUploaded(true);

      const filename = data.filename || data.split('/').pop();

      const downloadPath = `${import.meta.env.VITE_API_BASE_URL}api/download/${filename}`;
      setDownloadUrl(downloadPath);
    } catch (error) {
      console.error('Ошибка при загрузке:', error);
      toast.error('Ошибка при загрузке файла. Пожалуйста, попробуйте еще раз.');
    } finally {
      setIsLoading(false);
    }
  };

  const handleDownload = () => {
    if (downloadUrl) {
      const a = document.createElement('a');
      a.href = downloadUrl;
      a.download = '';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    }
  };

  return (
    <>
      <ToastContainer />
      <div className="upload-form">
        <input type="file" onChange={handleFileChange} />
        <button onClick={handleUpload} disabled={isLoading}>
          {isLoading ? 'Загрузка...' : 'Загрузить'}
        </button>
      </div>
      {isUploaded && (
        <div className="upload-form download">
          <p>Файл успешно загружен и обработан!</p>
          {downloadUrl && (
            <button onClick={handleDownload}>
              Скачать результат
            </button>
          )}
        </div>
      )}
    </>
  );
};

export default UploadForm;
