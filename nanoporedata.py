# -*- coding: utf-8 -*-
"""

@author: Mathieu Gibeault - Coda Consulting
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os
import ntpath
from glob import glob
import ConfigParser
import collections



class NanoporeData:
      
    def __init__(self):
        
        self.data = {}
        self.units = {}
        self.filters = {} 
        plt.style.use('default')
        mpl.rcParams['figure.figsize'] = [10.4, 7.15]
        mpl.rcParams['axes.labelsize'] = 14.3
        mpl.rcParams['axes.titlesize'] = 15.6
        mpl.rcParams['xtick.labelsize'] = 13
        mpl.rcParams['ytick.labelsize'] = 13            
        mpl.rcParams['legend.fontsize'] = 13        
          
    def load_data(self,DataPath, include=[], exclude=[], additionalParm = ''):
        
        files = getfilelist(DataPath,'PSD_Results*.csv',include,exclude)
        if not files:
            print('No PSD results file found in: ' + DataPath)
        else:
            self.data['PSD'] = {}
            for filepath in files:
                chip = getchipID(filepath)
                print("Reading file: " + filepath)
                df = pd.read_csv(filepath,skiprows=[1])
                #Removing null rows  (just checking one field, Date) and row where Bad_Data is not 0
                df = df[df['Date'].notnull() & (df['Bad_Data'] == 0)]          
                if not df.empty:
                    self.data['PSD'][chip] = df
            self.units['PSD'] = pd.read_csv(filepath,nrows=1,dtype=str).transpose()[0]
            self.units['PSD'][pd.isnull(self.units['PSD'])] = ''
        
        files = getfilelist(DataPath,'IV_Results*.csv',include,exclude)
        if not files:
            print('No IV results file found in: ' + DataPath)
        else:
            self.data['IV'] = {}
            for filepath in files:
                chip = getchipID(filepath)
                print('Loading ' + ntpath.basename(filepath) + " ...")
                df = pd.read_csv(filepath,skiprows=[1])
                #Removing null rows (just checking one field, Date) and rows where Pore size = NaN and where Bad_Data is not 0
                df = df[df['Date'].notnull() & df['Pore_Size_np'].notnull() &
                        df['Pore_Size_n'].notnull() & df['Pore_Size_p'].notnull() &
                          (df['Bad_Data'] == 0)]
                if not df.empty:
                    self.data['IV'][chip] = df
            self.units['IV'] = pd.read_csv(filepath,nrows=1,dtype=str).transpose()[0]
            self.units['IV'][pd.isnull(self.units['IV'])] = ''
        
        files = getfilelist(DataPath,'Fab_Stats*.ini',include,exclude)
        if not files:
            print('No Fab Stats file found in: ' + DataPath)
        else:
            self.data['fab_stats'] = {}
            for filepath in files:
                config = ConfigParser.RawConfigParser()
                config.optionxform = str #this is to make sure it does not change the case of option names
                config.read(filepath)
                self.data['fab_stats'][getchipID(filepath)] = ini_to_df('Values', config)
            self.units['fab_stats'] = ini_to_df('Units', config).transpose()[0]
        
            setrelativetime(self.data['PSD'],self.data['fab_stats'])
            setrelativetime(self.data['IV'],self.data['fab_stats'])
            
            if additionalParm != '':
                try:
                    df = pd.read_csv(additionalParm, index_col = 'PoreID')
                    for chip, fab_data in self.data['fab_stats'].items():
                        if chip in df.index:
                            for col in df:
                                self.data['fab_stats'][chip].insert(0,col,df.loc[chip][col])
                        else:
                            print(chip + ' Not in ' + additionalParm)
                except:
                    print("Could not read Additionnal Parameter file: " + additionalParm)
            
            self.exclude_chips()
            
            self.fab_stats_calculations()
            self.merge_PSD_IV()
            self.merge_fab_stats_data()
        
    
    def fab_stats_calculations(self):
        for chip, fab_data in self.data['fab_stats'].items():
            if 'Breakdown_Time' in fab_data:
                fab_data['Breakdown_Time'] = pd.to_datetime(fab_data['Breakdown_Time'][0], dayfirst=True)
            if ('Final_Pore_Size' in fab_data) and ('Post_Fab_Pore_Size' in fab_data):
                fab_data['Pore_Growth'] = fab_data['Final_Pore_Size'][0] - fab_data['Post_Fab_Pore_Size'][0]
                if 'Target_Size' in fab_data:
                    fab_data['Final_to_Desired_Pore_Size_Diff'] = fab_data['Final_Pore_Size'][0] - fab_data['Target_Size'][0]
                    fab_data['Size_Error'] = fab_data['Final_to_Desired_Pore_Size_Diff'][0] / fab_data['Target_Size'][0] * 100
                else: #this else statment was added because there were no Target Size value in earlier files and the target size was usually 5.
                    fab_data['Final_to_Desired_Pore_Size_Diff'] = fab_data['Final_Pore_Size'][0] - 5
                    fab_data['Size_Error'] = fab_data['Final_to_Desired_Pore_Size_Diff'][0] / 5 * 100
                if 'Cumulative_Energy' in fab_data:
                    fab_data['Cum_Energy_per_Growth'] = fab_data['Cumulative_Energy'][0] * 10**6 / fab_data['Pore_Growth'][0]                
                else:
                    print('Cumulative_Energy not found in fab_stats for chip# ' +  chip)
            else:
                print('Post_Fab or Final Pore Size not found in fab_stats for chip# ' +  chip)
            self.units['fab_stats']['Pore_Growth'] = 'nm'
            self.units['fab_stats']['Size_Error'] = '%'
            self.units['fab_stats']['Cum_Energy_per_Growth'] = 'uJ/nm'
            self.units['fab_stats']['Final_to_Desired_Pore_Size_Diff'] = 'nm'

    
    def exclude_chips(self):
        exclude_list = []
        for chip, fab_data in self.data['fab_stats'].items():
            if ('Exclude' in fab_data) and (fab_data['Exclude'][0] == 1): 
                exclude_list.append(chip)
        print('Excluding chips: ' + str(exclude_list))
        for dataset in ('PSD','IV','fab_stats'):
            for chip in exclude_list:
                if chip in self.data[dataset]:
                    del self.data[dataset][chip]
                
       
    def merge_fab_stats_data(self):
        for dataset in ('PSD','IV','PSD_IV'):
            for chip, fab_data in self.data['fab_stats'].items():
                if chip in self.data[dataset]:
                    for key in fab_data:
                        try:
                            self.data[dataset][chip].insert(0,key,fab_data[key][0])
                        except: 
                            pass
            self.units[dataset] = self.units[dataset].append(self.units['fab_stats'])
                    
    
    def merge_PSD_IV(self):
        self.data['PSD_IV'] = {}
        if (len(self.data['PSD']) <= len(self.data['IV'])):
            shortest_dict = self.data['PSD']
        else:
            shortest_dict = self.data['IV']
        for chip, data in shortest_dict.items():
            df = pd.merge(self.data['PSD'][chip],self.data['IV'][chip], 
                       how='inner', on='Iteration', suffixes=('', '_IV'))
            if not df.empty:
                self.data['PSD_IV'][chip] = df
                max_delta_t = (self.data['PSD_IV'][chip]['Relative_Time'] - self.data['PSD_IV'][chip]['Relative_Time_IV']).abs().max()
                if max_delta_t > 10:
                    print('Warning: Relative time difference of {:.0f} min found in PSD-IV merged data for chip# {}'.format(max_delta_t,chip))
            else:
                print('No PSD_IV merged data for chip: ' + chip)
                    
                    
        self.units['PSD_IV'] = self.units['PSD'].append(self.units['IV'])
        duplicates = self.units['PSD_IV'].index.get_duplicates()
        for key in duplicates:
            self.units['PSD_IV'][key+'_IV'] = self.units['PSD_IV'][key][1]
        self.units['PSD_IV'] = self.units['PSD_IV'].groupby(self.units['PSD_IV'].index).first()

    

    def make_chip_groups(self,group_by,dataset, name_in_label = True):
        groups = {}
        for chip, fab_data in self.data['fab_stats'].items(): 
            if group_by in fab_data:
                value = fab_data[group_by][0]
                if isinstance(value, str) and not value.isdigit():
                    label = value
                else:
                    if np.isnan(value):
                        continue
                    label = '{0:g}'.format(value)
                if label in groups:
                    groups[label].append(chip)
                else:
                    groups[label] = [chip]
        sorted_groups = collections.OrderedDict({})
        for key in sorted(groups.iterkeys()):
            if name_in_label:
                new_key = group_by + ' ' + key
            else:
                new_key = key
            sorted_groups[new_key] = groups[key]        
        #create new groups of data
        for type_group in sorted_groups:
            if not(self.filter_by_chipID(sorted_groups[type_group],dataset,type_group)):
                del sorted_groups[type_group]
        return sorted_groups
        
    def filter_data(self,query,data_name,filtered_data_name, how='all'):
        self.data[filtered_data_name] = {}
        self.units[filtered_data_name] = self.units[data_name]
        self.filters[filtered_data_name] = query
        for chip, data in self.data[data_name].items():
            try:
                if how =='first':
                    df = data.drop(data.index[1:]).query(query)
                else:
                    df = data.query(query)
                if not df.empty:
                    self.data[filtered_data_name][chip] = df
            except:
                print('Warning: Could not apply filter ( ' + query + ' ) to ' + chip)
            
    ### filter_by_fabstat not used anymore since fabstat data can be found in IV and PSD results datasets
    def filter_by_fabstat(self,query,data_name,filtered_data_name):
        self.data[filtered_data_name] = {}
        self.units[filtered_data_name] = self.units[data_name]
        self.filters[filtered_data_name] = query
        for chip, data in self.data['fab_stats'].items():
            try:
                if not data.query(query).empty:
                    if chip in self.data[data_name]: 
                        self.data[filtered_data_name][chip] = self.data[data_name][chip]
            except:
                print('Warning: Could not apply filter ( ' + query + ' ) to ' + chip)
                
    def filter_by_chipID(self,chipIDs,dataset_name,filtered_dataset_name):
        d = dict((i,self.data[dataset_name][i]) for i in self.data[dataset_name] if i in chipIDs)
        if d:
            self.data[filtered_dataset_name] = d
            self.units[filtered_dataset_name] = self.units[dataset_name]
            return True
        else:
            return False
            
    def plot_data(self,fields,data_name, grid=True, **kwargs):
        
        # the following code allows the 'fields' argument to be a string or a list
        if not isinstance(fields, (list, tuple)):
            fields = [fields]
            
        #setup the color scheme for the legend to make sure colors don't repeat
        NUM_COLORS = len(self.data[data_name])
        cm = plt.get_cmap('nipy_spectral')
        cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
        scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

        for Y in fields:  
            fig, ax = plt.subplots()
            ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
            X = 'Relative_Time'
            plt.xlabel('time [Min]')
            for chip, data in self.data[data_name].items():
                ax.plot(data[X], data[Y], label= chip, marker='o',**kwargs)
            ax.legend(loc=0, ncol=2)    
            if (Y in self.units[data_name]) and self.units[data_name][Y]:
                units = ' [' + self.units[data_name][Y] + ']'
            else:
                units = ''          
            plt.ylabel(Y.replace("_", " ") + units)
            if data_name in self.filters:
                plt.title(Y.replace("_", " ") + ' versus ' + X.replace("_", " ") + '\n' + 'Filter:  ' + self.filters[data_name])
            else:
                plt.title(Y.replace("_", " ") + ' versus ' + X.replace("_", " "))
            plt.legend(frameon=True)
            if grid:
                ax.grid(color='dimgrey',linewidth=0.5)
            plt.show() 
    
    def plot_hist(self,fields, data_sets, how = 'median', labels = '', bins = 100, 
                  log = False, freq=False, freq_overall=False, grid=True,  **kwargs):
        
        # the following code allows the 'fields' argument to be a string or a list
        if not isinstance(fields, (list, tuple)):
            fields = [fields]
        if not isinstance(data_sets, (list, tuple)):
            data_sets = [data_sets]
            
        for field in fields:  
            fig, ax = plt.subplots()
            hist_data_sets = []
            weight_sets = []
            total_len = 0
            for data_name in data_sets:
                hist_data = np.array([])
                for chip, data in self.data[data_name].items():
                    if field in data:
                        if how == 'all':
                            hist_data = np.append(hist_data, data[field].values)
                        elif how == 'first':
                            hist_data = np.append(hist_data, data[field].iloc[0])
                        elif how == 'last':
                            hist_data = np.append(hist_data, data[field].iloc[-1])
                        elif how == 'mean':
                            hist_data = np.append(hist_data, np.mean(data[field].values))
                        elif how == 'median':
                            hist_data = np.append(hist_data, np.median(data[field].values)) 
                if (hist_data.size != 0):
                    hist_data_sets.append(hist_data)
                    weight_sets.append(100 * np.ones_like(hist_data.astype(float)) / len(hist_data))
                    total_len = total_len + len(hist_data)
                else:
                    print('plot_hist -- No data for dataset: ' + data_name)
            if freq_overall:
                weight_sets[:] = [weight * len(data) / total_len for data, weight in zip(hist_data_sets,weight_sets)]            
            if labels == '':
                labels = data_sets
            if log == True:
                min_x = min(min(l) for l in hist_data_sets)
                max_x = max(max(l) for l in hist_data_sets)
                bins= np.logspace(np.log10(min_x),np.log10(max_x),bins)
                ax.set_xscale("log")
            if freq == True:
                ax.hist(hist_data_sets,label=labels, bins=bins, weights=weight_sets, **kwargs)
                plt.ylabel('frequency')
            else:
                ax.hist(hist_data_sets, label=labels, bins=bins, **kwargs)
                plt.ylabel('count')
            ax.legend(loc=0, ncol=2)    
            if (field in self.units[data_name]) and self.units[data_name][field]:
                units = ' [' + self.units[data_name][field] + ']'
            else:
                units = ''          
            plt.xlabel(field.replace("_", " ") + units)
            if (len(data_sets)==1) and (data_name in self.filters):
                plt.title(' Histogram of ' + field.replace("_", " ") + '\n' + 'Filter:  ' + self.filters[data_name])
            else:
                plt.title(' Histogram of ' + field.replace("_", " "))
            plt.legend(frameon=True)
            if grid:
                ax.grid(color='dimgrey',linewidth=0.5)
            plt.show()      
    
    
    def plot_2dhist(self,X,Y,data_set, how='all', bins=100, scale='linear', grid=True,
                xlim='', ylim=''):

        if not isinstance(how, (list, tuple)):
            how = [how,how]
        fig, ax = plt.subplots()
        X_data = np.array([])
        Y_data = np.array([])
        for chip, data in self.data[data_set].items():
            if (X in data) and (Y in data):
                if how[0] == 'all':
                    X_data = np.append(X_data, data[X].values)
                elif how[0] == 'first':
                    X_data = np.append(X_data, data[X].iloc[0])
                elif how[0] == 'last':
                    X_data = np.append(X_data, data[X].iloc[-1])
                elif how[0] == 'mean':
                    X_data = np.append(X_data, np.mean(data[X].values))
                elif how[0] == 'median':
                    X_data = np.append(X_data, np.median(data[X].values))
                    
                if how[1] == 'all':
                    Y_data = np.append(Y_data, data[Y].values)
                elif how[1] == 'first':
                    Y_data = np.append(Y_data, data[Y].iloc[0])
                elif how[1] == 'last':
                    Y_data = np.append(Y_data, data[Y].iloc[-1])
                elif how[1] == 'mean':
                    Y_data = np.append(Y_data, np.mean(data[Y].values))
                elif how[1] == 'median':
                    Y_data = np.append(Y_data, np.median(data[Y].values))
                           
        if not ylim == '':
            plt.ylim(ylim)   
        if not xlim == '':
            plt.xlim(xlim)
        ax.legend(loc=0, ncol=2)    
        if scale == 'logx':
            plt.xscale('log')
        elif scale == 'logy':
            bins= (bins , np.logspace(np.log10(min(Y_data)),np.log10(max(Y_data)),bins))
            plt.yscale("log")
        elif scale == 'loglog':
            plt.xscale('log')
            plt.yscale('log')
        
        plt.hist2d(X_data, Y_data,bins=bins)
        plt.colorbar() 
        if (Y in self.units[data_set]) and self.units[data_set][Y]:
            units = ' [' + self.units[data_set][Y] + ']'
        else:
            units = ''          
        plt.ylabel(Y.replace("_", " ") + units)
        if (X in self.units[data_set]) and self.units[data_set][X]:
            units = ' [' + self.units[data_set][X] + ']'
        else:
            units = ''          
        plt.xlabel(X.replace("_", " ") + units)
        plt.title(Y.replace("_", " ") + ' versus ' + X.replace("_", " "))
        if grid:
            ax.grid(color='dimgrey',linewidth=0.5)
        plt.show()    
    
    def plot_XY(self,X,Y,data_sets, how='all', labels='', scale='linear', grid=True,
                xlim='', ylim=''):

        if not isinstance(how, (list, tuple)):
            how = [how,how]
        if not isinstance(data_sets, (list, tuple)):
            data_sets = [data_sets]
        if labels == '':
            labels = data_sets 
        markers = ('.','d','*','p','o','<','H','P','X')    
        fig, ax = plt.subplots()
        
        for data_name, label, marker in zip(data_sets,labels,markers):
            X_data = np.array([])
            Y_data = np.array([])
            for chip, data in self.data[data_name].items():
                if (X in data) and (Y in data):
                    if how[0] == 'all':
                        X_data = np.append(X_data, data[X].values)
                    elif how[0] == 'first':
                        X_data = np.append(X_data, data[X].iloc[0])
                    elif how[0] == 'last':
                        X_data = np.append(X_data, data[X].iloc[-1])
                    elif how[0] == 'mean':
                        X_data = np.append(X_data, np.mean(data[X].values))
                    elif how[0] == 'median':
                        X_data = np.append(X_data, np.median(data[X].values))
                        
                    if how[1] == 'all':
                        Y_data = np.append(Y_data, data[Y].values)
                    elif how[1] == 'first':
                        Y_data = np.append(Y_data, data[Y].iloc[0])
                    elif how[1] == 'last':
                        Y_data = np.append(Y_data, data[Y].iloc[-1])
                    elif how[1] == 'mean':
                        Y_data = np.append(Y_data, np.mean(data[Y].values))
                    elif how[1] == 'median':
                        Y_data = np.append(Y_data, np.median(data[Y].values))
                        
            ax.plot(X_data, Y_data, label = label,marker=marker,linestyle='none')
        if not ylim == '':
            plt.ylim(ylim)   
        if not xlim == '':
            plt.xlim(xlim)
        ax.legend(loc=0, ncol=2)    
        if scale == 'logx':
            plt.xscale('log')
        elif scale == 'logy':
            plt.yscale('log')
        elif scale == 'loglog':
            plt.xscale('log')
            plt.yscale('log')
#        else:
#            plt.xscale('linear')
#            plt.yscale('linear')
        if (Y in self.units[data_name]) and self.units[data_name][Y]:
            units = ' [' + self.units[data_name][Y] + ']'
        else:
            units = ''          
        plt.ylabel(Y.replace("_", " ") + units)
        if (X in self.units[data_name]) and self.units[data_name][X]:
            units = ' [' + self.units[data_name][X] + ']'
        else:
            units = ''          
        plt.xlabel(X.replace("_", " ") + units)
        plt.title(Y.replace("_", " ") + ' versus ' + X.replace("_", " "))
        plt.legend(frameon=True)
        if grid:
            ax.grid(color='dimgrey',linewidth=0.5)
        plt.show()
        
    def plot_group(self,field, data_sets, field2 = '', how = 'median', scale = 'linear', grid=True,
                   ylim='', scale2='linear', y2lim='', how2 = 'median', **kwargs):
        
        fig, ax = plt.subplots() 
        if field2 != '':
            ax2 = ax.twinx()         
        labels = []
        idx = 0;
        for data_name in data_sets:
            Y_data = [] 
            Y2_data = []
            for chip, data in self.data[data_name].items():
                if (field in data):
                    if how == 'all':
                        Y_data.extend(data[field].values)
                    elif how == 'first':
                        Y_data.extend([data[field].iloc[0]])
                    elif how == 'last':
                        Y_data.extend([data[field].iloc[-1]])
                    elif how == 'mean':
                        Y_data.extend([np.mean(data[field].values)])
                    elif how == 'median':
                        Y_data.extend([np.median(data[field].values)])
                if (field2 != '') and (field2 in data):
                    if how2 == 'all':
                        Y2_data.extend(data[field2].values)
                    elif how2 == 'first':
                        Y2_data.extend([data[field2].iloc[0]])
                    elif how2 == 'last':
                        Y2_data.extend([data[field2].iloc[-1]])
                    elif how2 == 'mean':
                        Y2_data.extend([np.mean(data[field2].values)])
                    elif how2 == 'median':
                        Y2_data.extend([np.median(data[field2].values)])
                        
            if Y_data or Y2_data:
                labels.append(data_name)
                X_data = [idx]*len(Y_data)
                if (scale == 'logy') or (scale == 'log'):
                   ax.semilogy(X_data, Y_data, marker='o', color = 'b', linestyle='none')
                else:
                   ax.plot(X_data, Y_data, marker='o', color = 'b', linestyle='none') 
                if Y2_data:
                    X_data = [idx]*len(Y2_data)
                    if (scale2 == 'logy') or (scale2 == 'log'):
                        ax2.semilogy(X_data, Y2_data, marker='d', color = 'tab:orange', linestyle='none') 
                    else:
                        ax2.plot(X_data, Y2_data, marker='d', color = 'tab:orange', linestyle='none')
                idx += 1
        if grid:
            ax.grid(color='dimgrey',linewidth=0.5)
        plt.xticks(np.arange(len(labels)),labels)
        plt.xlim([-1,len(labels)])
        
        if ylim != '':
            ax.set_ylim(ylim)
        if (field in self.units[data_name]) and self.units[data_name][field]:
            units = ' [' + self.units[data_name][field] + ']'
        else:
            units = ''
        ax.set_ylabel(field.replace("_", " ") + units)                       
        if field2 != '': 
            if y2lim != '':
                ax2.set_ylim(y2lim)
            if (field2 in self.units[data_name]) and self.units[data_name][field2]:
                units = ' [' + self.units[data_name][field2] + ']'
            else:
                units = ''
            ax2.set_ylabel(field2.replace("_", " ") + units, color = 'tab:orange')  
            ax2.tick_params('y', colors='tab:orange', which='both')            
        plt.show()

      

    def split_PSD_by_bias(self,biases,delete_original):
        for chip, data in self.data['PSD'].items():
            for bias in biases:
                self.data['PSD'][chip + ' ' + str(bias) + 'mV'] = data[(data.Avg_Voltage < (float(bias)/1000)+0.05) & (data.Avg_Voltage > (float(bias)/1000)-0.05)]
            if delete_original:
                del self.data['PSD'][chip]



###############################################################################
#
#  Module Functions
#
###############################################################################

def getfilelist(path, searchString, include, exclude):  
    files = [y for x in os.walk(path) for y in glob(os.path.join(x[0], searchString))]
    results = []
    if include:
        for chip in include:
            results.extend([s for s in files if chip in s])
    else:
        results = files
    if exclude:
        for chip in exclude:
            results = [ x for x in results if chip not in x]
    return results

def getchipID(path):
    filenamestrings = ntpath.basename(path).split('_')
    return filenamestrings[2]

def setrelativetime(dataframes,fab_stats):
    for chip, data in dataframes.items():
        data['DateTime'] = pd.to_datetime(data['Date'].str.replace('-','/') + ' ' + data['Time'], dayfirst=True)
        try:
            StartTime = pd.to_datetime(fab_stats[chip]['Breakdown_Time'], dayfirst=True).iloc[0]
        except KeyError:
            print('No breakdown time for ' + chip + ' could be found.' )
            StartTime = data['DateTime'].iloc[0]
        data['Relative_Time'] = (data['DateTime'] - StartTime)/np.timedelta64(1, 'm')

def ini_to_dict(section, config):
    dictionary = {}
    for option in config.options(section):
        dictionary[option] = str_or_float(config.get(section, option).strip('"'))
    return dictionary

def ini_to_df(section, config):
    df = pd.DataFrame()
    for option in config.options(section):
        df[option] = [str_or_float(config.get(section, option).strip('"'))]
    return df

def str_or_float(s):
    try:
        return float(s)
    except ValueError:
        return s
     
def replaceNaN(unitseries):
    unitseries[unitseries == 0] = '' 

    


