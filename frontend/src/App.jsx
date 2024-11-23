import UploadForm from './components/UploadForm';
import ServiceDescription from './components/ServiceDescription';
import FutureFeatures from './components/FutureFeatures';
import './App.css';

const App = () => (
  <div className="app">
    <h1>Загрузка файла</h1>
    <UploadForm />
    <ServiceDescription />
    <FutureFeatures/>
  </div>
);

export default App;
