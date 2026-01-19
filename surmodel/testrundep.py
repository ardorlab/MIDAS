import multiprocessing as mp
import os
import pickle
import numpy as np
import time as TT
from temp import get_result, init_worker
# os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

if __name__ == "__main__":
    st = TT.time()
    # 1. Load shared data once in main process
    xsdict = pickle.load(open('/home/khnguy22/Deeponet-midas/MIDAS/surmodel/xsdata/axialFAtypeXS.pkl','rb'))
    fabulist = np.array([0, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.5, 15.0, 
                         17.5, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 
                         60.0, 65.0, 70.0, 75.0, 80.0])
    faaxial = np.array([15.24, 10.16, 5.08, 30.48, 30.48, 30.48, 30.48, 30.48,
                        30.48, 30.48, 30.48, 30.48, 30.48, 5.08, 10.16, 15.24])
    total_height = np.sum(faaxial)
    corebulist = [0.1, 0.4, 0.5, 1.0, 1.0] + [1.0]*28
    
    # 2. Define list of Loading Patterns (LPs)
    LPs = ['526  526  526  566  501  526  526  501  10 \n', 
       '526  526  462  502  526  566  501  501  10 \n', 
       '462  526  566  501  526  502  526  501  10 \n', 
       '566  462  526  526  501  526  462  501  10 \n', 
       '526  501  566  526  501  526  566  10   10 \n',
       '462  566  526  502  462  526  501  10   00 \n',
       '526  462  501  462  526  526  10   10   00 \n',
       '526  526  501  501  10   10   10   00   00 \n',
       '10   10   10   10   10   00   00   00   00 \n']

    LPs_o = " ".join(LPs).strip().split()
    list_of_LPs = [ LPs_o] 
    
    # 3. Setup Pool with initialization function
    with mp.Pool(processes=1, initializer=init_worker) as pool:
        # Prepare arguments
        args = [(lp, corebulist, xsdict, fabulist, faaxial, total_height, [0], [1,2,3,4,5,6,7,9,18,27,36,40,54,63]) 
                for lp in list_of_LPs]
        
        results = pool.starmap(get_result, args)

    print("Parallel processing complete:", results)
    print(f'time cost {TT.time()-st:.3f}')



# from test import get_result
# # # import Pool 

# LPs = ['526  526  526  566  501  526  526  501  10 \n', 
#        '526  526  462  502  526  566  501  501  10 \n', 
#        '462  526  566  501  526  502  526  501  10 \n', 
#        '566  462  526  526  501  526  462  501  10 \n', 
#        '526  501  566  526  501  526  566  10   10 \n',
#        '462  566  526  502  462  526  501  10   00 \n',
#        '526  462  501  462  526  526  10   10   00 \n',
#        '526  526  501  501  10   10   10   00   00 \n',
#        '10   10   10   10   10   00   00   00   00 \n']

# LPs_o = " ".join(LPs).strip().split()
# print(LPs)
# LPs = ["526","526","526" ,"566","501","526","526", "10", "10",
# "526","526","462" ,"502","526","566","501", "10", "10",
# "462","526","566" ,"501","526","502","526", "10", "10",
# "566","462","526" ,"526","501","526","462", "10", "10",
# "526","501","566" ,"526","501","526","566", "10", "10",
# "462","566","526" ,"502","462","526","501", "10", "00",
# "526","462","501" ,"462","526","526", "10", "10", "00",
# "526","462","462" ,"502", "10", "10", "10", "00", "00",
#  "10", "10", "10" , "10", "10", "00", "00", "00", "00",]

# for i in range(81):
#     print(LPs[i], LPs_o[i])


# stop
# print(get_result(LPs))