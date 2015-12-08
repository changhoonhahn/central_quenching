'''

Plotting class objects

'''
import matplotlib.pyplot as plt
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 

class Plots(object): 
    
    def __init__(self, **kwargs): 
        '''
        Class that generally encompasses all of the plotting I will do for the CenQue project
        '''
        self.kwargs = kwargs

        prettyplot()
        self.pretty_colors = prettycolors()

        self.fig = plt.figure(figsize=[10, 10])
        self.sub = self.fig.add_subplot(1,1,1)

    def save_fig(self, file_name): 
        '''
        save figure  to file_name
        '''

        self.fig.savefig(file_name, bbox_inches = 'tight')

        return None
