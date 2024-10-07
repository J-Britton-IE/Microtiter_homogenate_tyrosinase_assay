# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 17:19:04 2021

@author: UCD
"""

import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
import activity_analysis as aa
import pandas as pd
# Note the matplot tk canvas import
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import xlsxwriter

# VARS CONSTS:
_VARS = {'window': False}


def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


def delete_figure(figure):
    figure.get_tk_widget().forget()
    plt.close()
# \\  -------- PYSIMPLEGUI -------- //


AppFont = 'Any 16'
sg.theme('LightGrey')

font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 5}

plt.rc('font', **font)

Parameters = [[sg.Frame('Parameters', [[sg.Text('Molar Ex-coef: '), sg.Input('', size = (5, 2), key = 'epsilon' ),
                                        sg.Text('Replicates: '), sg.Input('', size = (5, 2), key = 'replicateNumber'),
                                        sg.Text('Stepwidth: '), sg.Input('', size = (5, 2), key = 'stepWidth'),
                                        sg.Text('R2 limit: '), sg.Input('', size = (5,2), key = 'r2_lim')]], border_width =  5)]]

Rows       = [[sg.Frame('Row Conc.', [[sg.Text('B: '), sg.Input('', size = (5,2), key = 'Row_B')],
          [sg.Text('C: '), sg.Input('', size = (5,2), key = 'Row_C')],
          [sg.Text('D: '), sg.Input('', size = (5,2), key = 'Row_D')],
          [sg.Text('E: '), sg.Input('', size = (5,2), key = 'Row_E')],
          [sg.Text('F: '), sg.Input('', size = (5,2), key = 'Row_F')],
          [sg.Text('G: '), sg.Input('', size = (5,2), key = 'Row_G')],
          [sg.Text('H: '), sg.Input('', size = (5,2), key = 'Row_H')]], border_width = 5)]]


tab1 = [[sg.Canvas(key ='figCanvas')]]
tab2 = [[sg.Canvas(key = 'figCanvas2')]]
tab3 = [[sg.Canvas(key = 'figCanvas3')]]
tab4 = [[sg.Canvas(key = 'figCanvas4')]]
tab5 = [[sg.Canvas(key = 'figCanvas5')]]

buttons = [[sg.Frame('', [
           [sg.Button('Clear', use_ttk_buttons=True), 
           sg.Button('Help', use_ttk_buttons= True), 
           sg.Button('Exit', use_ttk_buttons= True)]], border_width = 5, element_justification='right')]]


data_outs  = [[sg.Frame('Data',
                        [[sg.TabGroup([[#sg.Tab('', tab0),
                                        sg.Tab('OD400', tab1), 
                                        sg.Tab('Rates', tab2), 
                                        sg.Tab('Linearity', tab3),
                                        sg.Tab('R2 values', tab4), 
                                        sg.Tab('Linear rates', tab5)
                  ]],)],
                [sg.Text('1ug/mM/min max rate (mM/min):'), sg.Input('', size = (15,2), key ='max_result')]
                ], border_width = 5)]]

analyse_plot = [[sg.Frame('', [[sg.Button('Analyze', use_ttk_buttons= True), sg.Button('Plot', use_ttk_buttons= True)]], border_width= 5)]]

savings   = [[sg.Frame('',[[sg.Text('Results file name: '), sg.Input('', size=(55, 2), key = 'resultfilename')],
          
          [sg.FolderBrowse('Save folder', size=(10, 1)), sg.Input(size=(58, 1), key='resultfolder')], 
          
          [sg.Button('Save')]], border_width = 5)]]


layout = [[sg.FileBrowse('Select File', size=(10, 1)), sg.Input('', size=(60, 2), key ='myFile')],
          
          [sg.Column(Parameters, element_justification='c')],
          
          [sg.Column(Rows, element_justification='left', vertical_alignment = 'top'), 
           sg.Column(data_outs, element_justification='left', vertical_alignment='top')],
          
          [sg.Column(analyse_plot, element_justification= 'left')],
          
          [sg.Column(savings, element_justification='left', vertical_alignment = 'top')],
          
          
          [sg.Column(buttons, element_justification='right', vertical_alignment= 'bottom')]]
         
                                                                    
_VARS['window'] = sg.Window('JBio - Colourometric enzyme activity analysis',
                            layout,
                            finalize=True,
                            resizable=True)
                            # ,icon=r'C:\Users\UCD\PycharmProjects\Bioinformatics\96_well_plate_2.png').read(close=True)
                            #size = (540, 550))


# MAIN LOOP
while True:
    event, values = _VARS['window'].read(timeout=200)
    
    results = []  
        
    if event == 'Analyze':
        
        if values['myFile'] == '':
            sg.Popup('Please ensure a datafile has been slected',
                 title='Warning - No datafile',
                 keep_on_top= True)
        
        elif values['epsilon'] == '':
            sg.Popup('Please enter Molar extinction coefficient.',
                     title ='Warning!',
                     keep_on_top= True)
        
        elif values['replicateNumber'] == '':
            sg.Popup('Please enter number of replicates.',
                     title ='Warning!',
                     keep_on_top= True)
            
        elif values['stepWidth'] == '':
            sg.Popup('Please select stepwidth for rate calculations.',
                     title ='Warning!',
                     keep_on_top= True)
        
        elif values['r2_lim'] == '':
            sg.Popup('Please enter an R2 limit from 0.9 to 1')
            
        elif values['Row_B'] == '' and values['Row_C'] == '' and values['Row_D'] == '' and values['Row_E'] == '' and values['Row_F'] == '' and values['Row_G'] == '' and values['Row_H'] == '':
                sg.Popup('Please enter protein loading per row.',
                         title = 'Warning!',
                         keep_on_top= True)
        
        elif values['myFile'] != '':
            
            myFile  = values['myFile'];
            myFile  = pd.read_excel(myFile, header = None)
        
            epsilon = int(values['epsilon']);
            
            r2_limit = float(values['r2_lim'])
            
            if r2_limit > 1:
                sg.Popup('R2 value above 1, please enter a correct value.')
                continue
        
            EnzymeConcentrations = [];
            EnzymeConcentrations.append(values['Row_B'])
            EnzymeConcentrations.append(values['Row_C'])
            EnzymeConcentrations.append(values['Row_D'])
            EnzymeConcentrations.append(values['Row_E'])
            EnzymeConcentrations.append(values['Row_F'])
            EnzymeConcentrations.append(values['Row_G'])
            EnzymeConcentrations.append(values['Row_H'])
            EnzymeConcentrations = [x for x in EnzymeConcentrations if x != '']
            EnzymeConcentrations = [float(i) for i in EnzymeConcentrations]
                
            replicates        = int(values['replicateNumber']);
        
            stepWidth         = int(values['stepWidth']);
            
            try:
                results           = aa.activity_analysis(epsilon, EnzymeConcentrations, replicates, stepWidth, myFile, r2_limit);
            except:
                sg.Popup('Unable to analyze.\n'
                         'Check input file and input parameters.',
                         title = 'Warning!',
                         keep_on_top= True)
                continue
                
            time              = results[0]
            mean_data         = results[1]
            rate_times        = results[2]
            reaction_rates    = results[3]
            linearity         = results[4]
            linearity         = linearity.iloc[:, 1:]
            rsq_plot_vals     = results[5]
            rsq_plot_vals     = rsq_plot_vals.iloc[:, 1:]
            max_lin_rates     = results[6]
            max_lin_rate_time = max_lin_rates.iloc[:, 0]
            max_lin_rates     = max_lin_rates.iloc[:, 1:]
            max_lin_rate_time = max_lin_rate_time.to_list()
            max_lin_rate_stdev = results[7]
            max_lin_rate_stdev = max_lin_rate_stdev.to_list()
            max_lin_rates       = max_lin_rates.values.tolist()
            max_lin_rates       = sum(max_lin_rates, [])
            EZC2 = EnzymeConcentrations
            EZC2.insert(0, 0)
            legends = EnzymeConcentrations[1:]
        
    elif event == 'Plot':
        
        
        fig = plt.figure()
        plt.plot(time, mean_data, linestyle ='solid')
        plt.ylabel('OD400')
        plt.xlabel('Time (mins)')
        plt.title('OD400 over time')
        plt.legend(legends)
        fig.set_size_inches(4,2.5)
        draw_figure(_VARS['window']['figCanvas'].TKCanvas, fig)
        # plt.savefig('OD400_plot.png')
        
        fig2 = plt.figure()
        plt.plot(rate_times, reaction_rates, linestyle ='solid')
        plt.ylabel('Product formation mM/min')
        plt.xlabel('Time (mins)')
        plt.title('Reaction rate mM/min')
        plt.legend(legends)
        fig2.set_size_inches(4,2.5)
        draw_figure(_VARS['window']['figCanvas2'].TKCanvas, fig2)
    
        fig3 = plt.figure()
        plt.plot(EZC2, linearity, linestyle ='solid')
        plt.title('Protein loading vs Rate at each timepoint')
        plt.ylabel('Product formation mM/min')
        plt.xlabel('Protein loading ug/mM/ml')
        fig3.set_size_inches(4,2.5)
        draw_figure(_VARS['window']['figCanvas3'].TKCanvas, fig3)
        
        fig4 = plt.figure()
        plt.plot(rate_times, rsq_plot_vals, linestyle ='solid')
        plt.ylabel('R2 value')
        plt.xlabel('Time (mins)')
        plt.title('R2 value over time')
        fig4.set_size_inches(4,2.5)
        draw_figure(_VARS['window']['figCanvas4'].TKCanvas, fig4)
    
        fig5 = plt.figure()
        y_pos = np.arange(len(max_lin_rate_time))
        
        plt.bar(y_pos, max_lin_rates, yerr = max_lin_rate_stdev, align = 'center', alpha = 0.5)
        plt.ylim((0, 0.1))
        plt.xticks(y_pos, max_lin_rate_time)
        plt.xlabel('Timepoints (mins)')
        plt.ylabel('Product formation mM/min')
        plt.title('Timepoints where protein loading vs Rate is linear.')
        fig5.set_size_inches(4,2.5)
        draw_figure(_VARS['window']['figCanvas5'].TKCanvas, fig5)
        
        
        max_res_str = max(max_lin_rates)
        
        max_res_str = str(round(max_res_str, 5)) + 'mM/min'
        
        _VARS['window']['max_result'](max_res_str)
        
    elif event == 'Help':
        sg.Popup('Microtiter plate ctivity analysis tool.\n'
                 'Ensure plate set up as follows: \n'
                 ' - Row A = Blanks ONLY. \n'
                 ' - Every row should contain 1 protein concentration. \n'
                 ' - Columns 1 - 12 should contain replicates. \n',
                 title='Help',

                 keep_on_top= True)
    
    elif event == 'Clear':
        _VARS['window']['myFile']('')
        _VARS['window']['epsilon']('')
        _VARS['window']['replicateNumber']('')
        _VARS['window']['stepWidth']('')
        _VARS['window']['Row_B']('')
        _VARS['window']['Row_C']('')
        _VARS['window']['Row_D']('')
        _VARS['window']['Row_E']('')
        _VARS['window']['Row_F']('')
        _VARS['window']['Row_G']('')
        _VARS['window']['Row_H']('')
        _VARS['window']['max_result']('')
        _VARS['window']['resultfilename']('')
        _VARS['window']['resultfolder']('')
        _VARS['window']['r2_lim']('')
         
        
    elif event == 'Save':
        
        results_location = values['Save folder']
    
        _VARS['window']['resultfolder'](results_location)
      
        save_folder = str(values['resultfilename'])
        
        save_file_path = str(results_location) + "/" + str(save_folder)
        
        save_spot = str(results_location) + "/" + str(save_folder) + "/" + str(values['resultfilename'] + '.xlsx')
        
        os.mkdir(save_file_path)
        
        resultsFileName = str(values['resultfilename']) + '.xlsx'
        
        writer = pd.ExcelWriter(save_spot, engine='xlsxwriter')
        
        myFile.to_excel(writer, sheet_name = 'Raw_data', header = False, index = False)
        
        mean_data_time = pd.concat([time, mean_data],axis = 1)
        
        mean_data_time.to_excel(writer, sheet_name = 'Mean OD400', header = True, index = False )
        
        rates_with_times = pd.concat([rate_times, reaction_rates], axis =1)
        
        rates_with_times.to_excel(writer, sheet_name = 'Rates', header = True, index = False)
        
        linearity_at_timepoints = pd.concat([rate_times, rsq_plot_vals], axis = 1)
        
        linearity_at_timepoints.to_excel(writer, sheet_name = 'R2 values', header = True, index = False)
        
        max_lin_rate_time = pd.DataFrame(max_lin_rate_time, columns = ['Time'] )
        
        max_lin_rates = pd.DataFrame(max_lin_rates, columns = ['Rate'])
        
        max_lin_rate_stdev = pd.DataFrame(max_lin_rate_stdev, columns = ['STDEV'])
        
        max_linear_rates_excel = pd.concat([max_lin_rate_time, max_lin_rates, max_lin_rate_stdev], axis = 1)
        
        max_linear_rates_excel.to_excel(writer, sheet_name = 'Linear rates', header = True, index = False)
        
        writer.save() 
       
        fig.savefig(str(save_file_path) + '/' + '1_OD400_plot.png')
        fig2.savefig(str(save_file_path) + '/' + '2_Rates_plot.png')
        fig3.savefig(str(save_file_path) + '/' + '3_linearity_plot.png')
        fig4.savefig(str(save_file_path) + '/' + '4_R2_plot.png')
        fig5.savefig(str(save_file_path) + '/' + '5_Max_linear_rates_plot.png')
        
    elif event == sg.WIN_CLOSED or event == 'Exit':
        break
    
_VARS['window'].close()
