import numpy as np 

def pyreplicate(dict, n): 
    '''
    Parameters
    ----------
    dict : 
        Python dictionary that specifies the type of data 
    n : 
        Number of times you want it replicated
    '''
    
    for key in dict.keys(): 
        print type(dict[key])
        dict[key] = np.repeat(dict[key], n)
    return dict

if __name__=='__main__': 
    dict = {'ra': 0.0, 'dec': 0.0, 'z': 0.0, 'blah': 'asdf'}

    print pyreplicate(dict, 2)
