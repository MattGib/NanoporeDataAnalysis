import pandas as pd
import numpy as np
from nanoporedata import NanoporeData
import matplotlib.pyplot as plt
import collections


include= [] # ex: include =  ['AMJ05','AMJ04','AMJ06','SK186']
exclude= [] #ex: exclude = ['AMJ016','AMJ017']


nd = NanoporeData()
# load_data loads all data from specified path.  
# Creates data['PSD'] data['IV'] and data['PSD_IV']
nd.load_data('C:\\Data\\Result Files',include=include,exclude=exclude, 
             additionalParm = 'C:\\Data\\Result Files\\Additional Parameters.csv' )


#create type groups and data sets
type_groups  = nd.make_chip_groups('Type','fab_stats',name_in_label=True)


nd.filter_data('(0.19<Avg_Voltage<0.21) and (L<=0.003)','PSD', 'Great',how='first')
nd.filter_data('(0.19<Avg_Voltage<0.21)  and (0.003<L<=0.01)','PSD', 'Good',how='first')
nd.filter_data('(0.19<Avg_Voltage<0.21)  and (L>0.01)','PSD', 'Poor',how='first')


nd.plot_hist('L',['Great','Good','Poor'],how='first',bins=3,log=False,freq=True,color=['maroon','blue','red'])
nd.plot_hist('L',['Great','Good','Poor'],how='first',bins=1,log=False,freq=True,freq_overall=True)

nd.plot_hist('L',['Great','Good','Poor'],how='first',bins=10,log=False,freq=False)
ax = plt.gca()
ax.set_yticks([0,10,20])
ax.set_yticks([0,5,10,15,20], minor = True)


#Fabrication Type vs Initial Pore Size
nd.plot_group('Post_Fab_Pore_Size',type_groups.keys(),how ='first',scale = 'linear')

#Fabrication Type vs Cumulative enery per nm of pore growth
nd.plot_group('Cum_Energy_per_Growth',type_groups.keys(),how ='first',scale = 'log',
              field2='Size_Ratio_pn', how2 = 'median', scale2 = 'linear')

#Fabrication Type vs Conditioning Time
nd.plot_group('Conditioning_Time',type_groups.keys(),how ='first',scale = 'linear',ylim=(0,80))


nd.plot_group('Final_Pore_Size',type_groups.keys(),how ='first',scale = 'linear')


#Fabrication Type vs Difference between Desired and Final pore size
nd.plot_group('Final_to_Desired_Pore_Size_Diff',type_groups.keys(),how='first',scale = 'linear')

nd.plot_XY('Final_to_Desired_Pore_Size_Diff','Conditioning_Time',type_groups.keys(),ylim=(0,80),grid=False)
nd.plot_XY('Size_Error','Conditioning_Time',type_groups.keys(),ylim=(0,80),grid=False)



nd.filter_data('(-0.19 > Avg_Voltage > -0.21) and (4 < Final_Pore_Size < 7)','PSD', 'PSD -200mV')
nd.filter_data('(0.19 < Avg_Voltage < 0.21) and (4 < Final_Pore_Size < 7)','PSD', 'PSD +200mV')
nd.filter_data('(4 < Final_Pore_Size < 7) and ((-0.19 > Avg_Voltage > -0.21) or (0.19 < Avg_Voltage < 0.21))','PSD', 'Size4to7')

nd.plot_XY('Breakdown_Time','L',['PSD -200mV','PSD +200mV'],how =['first','first'],scale='logy')

nd.plot_XY('Baseline_Current','L',['PSD -200mV','PSD +200mV'],how =['first','first'])

nd.plot_2dhist('Baseline_Current','L','PSD +200mV',how =['first','first'],bins=15,scale = 'logy')
nd.plot_2dhist('Baseline_Current','L','Size4to7',how =['first','first'],bins=15)


del nd.data['Type 1']['AMJ09'] # remove because conditioning time is wrong
nd.plot_hist('Conditioning_Time', ('Type 1','Type 8'),bins=50, freq=False)

R_groups  = nd.make_chip_groups('R','PSD +200mV',name_in_label=True)
nd.plot_hist('L', R_groups.keys(),bins=20, log=True, freq=True)
nd.plot_XY('Final_Pore_Size','L',R_groups.keys(),how =['first','median'],scale='logy')


#nd.plot_hist('L', ('Delay -200mV','No Delay -200mV'),bins=12, log=True, freq=True)
#nd.plot_hist('L', ('Delay -200mV','No Delay -200mV'),bins=12, log=True, freq=True, histtype='stepfilled')
#nd.plot_hist('L', ('Delay -200mV','No Delay -200mV'),bins=12, log=True, freq=True, histtype='step', linewidth=2)


