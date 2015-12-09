'''

SMF plotting class 

'''
from smf import SMF
from plots import Plots
from treepm.sham import SMFClass 

class PlotSMF(Plots): 

    def __init__(self, cenque=None, **kwargs): 
        ''' 
        Child class of Plots class that describes SMF plots for 
        different class objects including CenQue and GroupCat
        '''
        super(PlotSMF, self).__init__(**kwargs)

        self.fig = plt.figure(figsize=(10,10))
        self.sub = self.fig.add_subplot(111)

    def cenque(self, cenque, type='total', **kwargs): 
        ''' 
        Plot SMF distribution for CenQue class object
        
        Parameters
        ----------
        cenque : 
            Central Quenching Class Object
        type : 
            type of SMF plot. If type = 'total', then plot the total SMF. 
            if tyep='sfq' then plot the SF and Q components separately in the
            SMF. 
        '''
        if 'label' in kwargs: 
            smf_label = kwargs['label']
            kwargs.pop('label', None)
        else: 
            smf_label = 'z ='+str(cenque.zsnap) 
    
        if 'line_color' in kwargs: 
            line_color = kwargs['line_color']
            kwargs.pop('line_color', None)
        else: 
            try: 
                line_color = self.pretty_colors[cenque.nsnap]
            except TypeError: 
                line_color = 'black'

        if 'line_style' in kwargs:
            line_style = kwargs['line_style'] 
            kwargs.pop('line_style', None)
        else: 
            line_style = '-'

        if 'lw' in kwargs: 
            line_width = kwargs['lw'] 
            kwargs.pop('lw', None)
        else:
            line_width = 4

        #mf = SMF()
        mass, phi = cenque.SMF(**kwargs) #mf.cenque(cenque, **mkwargs)

        self.sub.plot(mass, phi, 
                c=line_color, ls=line_style, lw=line_width,
                label=smf_label
                )
    
    def analytic(self, redshift): 
        '''
        Analytic SMF at redshift 
        '''
        if 'label' in kwargs: 
            smf_label = kwargs['label']
            kwargs.pop('label', None)
        else: 
            smf_label = 'Analytic'
    
        if 'line_color' in kwargs: 
            line_color = kwargs['line_color']
            kwargs.pop('line_color', None)
        else: 
            line_color = 'black'

        if 'line_style' in kwargs:
            line_style = kwargs['line_style'] 
            kwargs.pop('line_style', None)
        else: 
            line_style = '--'

        if 'lw' in kwargs: 
            line_width = kwargs['lw'] 
            kwargs.pop('lw', None)
        else:
            line_width = 3
        smf = SMF()
        mass, phi = smf.analytic(redshift) 
        
        self.sub.plot(mass, phi, 
                c=line_color, ls=line_style, lw=line_width,
                label=smf_label
                )

    def groupcat(self, Mrcut=18, **kwargs): 
        ''' 
        Plot sSFR distribution for Group Catalog data
        '''

        raise NotImplementedError

    def set_axes(self): 
        ''' 
        Set up axes of figures
        '''

        self.sub.set_ylim([10**-5, 10**-1])
        self.sub.set_xlim([8.75, 12.0])

        # x and y labels
        self.subs.set_xlabel(r'Mass $\mathtt{M_*}$') 
        self.subs.set_ylabel(r'Stellar Mass Function $\mathtt{\Phi}$', fontsize=20) 
        
        self.subs.legend(loc='upper left', frameon=False)

        return None

    def save_fig(self, file_name): 
        ''' 
        save figure to file 
        '''
        self.fig.savefig(file_name, bbox_inches='tight')
