import pickle as pkl
import numpy as np
import joblib

def convert(pickle_file_path):
    joblib_path = pickle_file_path.replace('.pkl', '.joblib')
    with open(pickle_file_path, 'rb') as f:
        data = pkl.load(f)
    if isinstance(data, (list, tuple, dict)):
        print('convertible???')
        for key in data:
            if isinstance(data[key], list):
                data[key] = np.array(data[key])
                print(f"Optimized key '{key}': converted list to numpy array.")
        
        # Save as joblib
        joblib.dump(data, joblib_path)
        print(f"Success! Dictionary saved to: {joblib_path}")
    else:
        print('can not??')

testpath = '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/base_model/traindataall_coredata/scaler.pkl'
convert(testpath)