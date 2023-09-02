class Data_Handler:
    def __init__(self, data):
        import numpy as np
        
        self.data = data
        
        global Data, n_passos, n_sitios
        
        Data = data
        Data = np.array(Data)
        
        n_passos = len(Data)
        n_sitios = int((len(Data[0])-1)/2)
        
        
    def handler(self):
        import numpy as np
        t = Data[:, 0]
        
        Gt1 = np.zeros((n_sitios, 2, n_passos))
        
        for i in range(0, n_sitios):
            Gt1[i][0] = Data[:, 1 + 2*i]
            Gt1[i][1] = Data[:, 2 + 2*i]
            
        
        return t, Gt1
    
    def prepross(self, lb, sequence, val_split = 0):
        import numpy as np
        
        X, y = list(), list()
        for i in range(len(sequence)):
            #ind the end of this pattern
            end_ix = i + lb
    		# check if we are beyond the sequence
            if end_ix > len(sequence)-1:
                break    		
            # gather input and output parts of the pattern
            seq_x, seq_y = sequence[i:end_ix], sequence[end_ix]
            X.append(seq_x)
            y.append(seq_y)
            
        X = np.array(X)
        y = np.array(y)
        
        X = X.reshape((len(X), lb, 1, 1))
        
        split_l = int(val_split*len(X))
        
        X_te = X[split_l:]
        y_te = y[split_l:]
        
        X_tr = X[:split_l]
        y_tr = y[:split_l]
        
        return (X_tr, y_tr), (X_te, y_te)
        
     
