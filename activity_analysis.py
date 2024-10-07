# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 15:02:16 2021

@author: UCD
"""
# microtitre activity assay analysis function
import pandas as pd
import split

def activity_analysis(epsilon, EnzymeConcentrations, replicateNumber, stepWidth, myFile, rsq_limit):


# Matlab data analysis script translation

    pathlength            = 0.622; # Based on cuvettes used
       
    myFile     = myFile.drop(myFile.index[0:10]); # Drop first 12 rows, not useful data.

    myFile.dropna(how = 'all', axis = 1, inplace = True); # Drop empty cells

    rawData    = myFile.drop(myFile.index[0:2]); # Numbers only

    rawData.drop(rawData.iloc[:,0:2], inplace = True, axis = 1); # Remove time column

    sampleData = rawData.drop(rawData.columns[0:3], axis =1); # Sample data only, no blanks

    sampleData = sampleData.drop(sampleData.index[0]);

    myGroups = myFile.values[1].tolist(); # make new group containing group ID (Row on plate)

    del myGroups[0:5];

    myGroupSet = set(myGroups);

    myGroupSet = list(myGroupSet);

    myGroupSet.sort();

    myGroupSet = ['Time'] + myGroupSet;

    time = myFile.drop(myFile.columns[0:1], axis = 1); # Get new array with time values.

    time = time.drop(time.columns[1::], axis = 1);

    time = time.drop(time.index[0:3]);

    time.columns = ['Times'];

    time_numeric = len(time);

    time_numeric_list = []

    for i in range(time_numeric):
        tee = (i * 0.5);
        time_numeric_list.append(tee)
    
    
    time_numeric_list = pd.DataFrame(time_numeric_list);

    time_numeric_list.columns = ['Times']

    time = time_numeric_list;

    blanks = rawData.drop(rawData.columns[3::], axis = 1); # Get blank values

    blanks = blanks.drop(blanks.index[0]);

    blank_means = blanks.mean(axis = 'columns');

    samples_blk_adj = sampleData.subtract(blank_means, axis = 'index'); # subtraxt blanks from sample data.

    productConcentrations = samples_blk_adj/(epsilon * pathlength); #convert OD value to uM

    productConcentrations_mM = productConcentrations * 1000; # convert to mM

    no_samples = len(productConcentrations_mM.columns);

    no_bio_samples = no_samples / replicateNumber;

    productConcentrations_mM.columns = myGroups;
    
    row_b = productConcentrations_mM.filter(regex = 'B');

    row_b_mean = row_b.mean(axis = 'columns');

    row_b_mean = pd.DataFrame(row_b_mean);

    row_b_mean.reset_index(inplace = True);

    row_b_mean = row_b_mean.iloc[:, 1:];

    row_c = productConcentrations_mM.filter(regex = 'C');

    row_c_mean = row_c.mean(axis = 'columns');

    row_c_mean = pd.DataFrame(row_c_mean);

    row_c_mean.reset_index(inplace = True);

    row_c_mean = row_c_mean.iloc[:, 1:];

    row_d = productConcentrations_mM.filter(regex = 'D');

    row_d_mean = row_d.mean(axis = 'columns');

    row_d_mean = pd.DataFrame(row_d_mean);

    row_d_mean.reset_index(inplace = True);

    row_d_mean = row_d_mean.iloc[:, 1:];

    row_e = productConcentrations_mM.filter(regex = 'E');

    row_e_mean = row_e.mean(axis = 'columns');

    row_e_mean = pd.DataFrame(row_e_mean);

    row_e_mean.reset_index(inplace = True);

    row_e_mean = row_e_mean.iloc[:, 1:];

    row_f = productConcentrations_mM.filter(regex = 'F');

    row_f_mean = row_f.mean(axis = 'columns');

    row_f_mean = pd.DataFrame(row_f_mean);

    row_f_mean.reset_index(inplace = True);

    row_f_mean = row_f_mean.iloc[:, 1:];

    row_g = productConcentrations_mM.filter(regex = 'G');

    row_g_mean = row_g.mean(axis = 'columns');

    row_g_mean = pd.DataFrame(row_g_mean);

    row_g_mean.reset_index(inplace = True);

    row_g_mean = row_g_mean.iloc[:, 1:];

    row_h = productConcentrations_mM.filter(regex = 'H');

    row_h_mean = row_h.mean(axis = 'columns');

    row_h_mean = pd.DataFrame(row_h_mean);

    row_h_mean.reset_index(inplace = True);

    row_h_mean = row_h_mean.iloc[:, 1:];


## Plotting of mean OD change per protein loading.

    mean_data = pd.concat([time, row_b_mean, row_c_mean, row_d_mean, row_e_mean, row_f_mean, row_g_mean, row_h_mean], axis = 1);

    mean_data.dropna(how = 'all', axis = 1, inplace = True);

    mean_data.columns = myGroupSet;

    OD_plot = mean_data.plot(x = 'Time');

    OD_plot.set_xlabel('Time (mins)');

    OD_plot.set_ylabel('OD400');

    OD_plot.legend(EnzymeConcentrations);

    OD_plot.set_title('OD400 change over time per protein loading (ug/mM/ml)');

    mean_data.drop('Time', inplace = True, axis = 1);

# Find reaction rates and plot

    time.columns = ['Time'];

    time_list = time['Time'].to_list();

    rates = split.rate_func(time_list, mean_data, stepWidth);

    rates.columns = myGroupSet;

    myZeros = [0] * int(len(myGroupSet))

    myZeros = pd.DataFrame(myZeros); 

    myZeros = myZeros.transpose();

    myZeros.columns = myGroupSet;

    rate_frames = [myZeros, rates];

    rate_zeros = pd.concat(rate_frames, axis = 0, ignore_index = True);

    rate_plot = rate_zeros.plot(x = 'Time');

    rate_plot.legend(EnzymeConcentrations);

    rate_plot.set_xlabel('Time (mins)');

    rate_plot.set_title('Reaction rates over time')

    rate_plot.set_ylabel('Product formation uM/min');

# Make linear graph of protein conc vs rate at each timepoint.

    transposed_rates = rate_zeros.transpose();

    transposed_rates.columns = transposed_rates.iloc[0];

    transposed_rates = transposed_rates.iloc[1: , :]

    transposed_rate_zeros = pd.DataFrame([[0]*transposed_rates.shape[1]], columns = transposed_rates.columns);

    transposed_rate_zeros.columns = transposed_rates.columns;

    transposing = [transposed_rate_zeros, transposed_rates]

    transposed = pd.concat(transposing, axis = 0, ignore_index = True);

    EnzymeConcentrations.insert(0,0);

    Enzyme_Cs = pd.DataFrame(EnzymeConcentrations, columns = ['Protein loading (ug/mM)']);
    
    linearity_concat = [Enzyme_Cs, transposed];

    linearity = pd.concat(linearity_concat, axis = 1, ignore_index = False);

    linearity_plot = linearity.plot(x = 'Protein loading (ug/mM)', legend = None);

    linearity_plot.set_ylabel('Product formation mM');

    linearity_plot.set_title('Linear correlation between protein loading and reaction rate at each time point');

# Next to do - r2 value at each timepoint and graph.
    r2_vals = []

#for i in range(0, len(linearity.columns)):
    
    for column in linearity.columns[1:]:
    
        column_list = linearity[column].tolist();
    
        score = split.rsquared(EnzymeConcentrations, column_list);
    
        r2_vals.append(score);
    
    rate_times = rate_zeros['Time'];

    rate_times = rate_times.to_frame()
    
    r2_vals    = pd.DataFrame(r2_vals, columns = ['R2 values']);

    r2_all_together = [rate_times, r2_vals]

    rsq_plot_vals = pd.concat(r2_all_together, axis = 1, ignore_index = False);

    rsq_plot = rsq_plot_vals.plot(x = 'Time');

    rsq_plot.set_xlabel('Time (mins)');

    rsq_plot.set_ylabel('R2 value');

    rsq_plot.set_title('Correlation between protein loading & Rate over time');

    #rsq_limit = 0.98;

    high_correlation_time = rsq_plot_vals[rsq_plot_vals['R2 values'] > rsq_limit];

    high_correlation_time = high_correlation_time['Time'];

    high_correlation_time = high_correlation_time.tolist();

# Have all the times at which rate to loading correlation is high. Now find out what those rates are and plot.
    high_correlation_time_rates = []

    rate_zeros = rate_zeros.set_index('Time')

    for index, rows in rate_zeros.iterrows():
        for item in high_correlation_time:
            if index == item:
                high_correlation_time_rates.append(rows)

    high_correlation_time_rates = pd.DataFrame(high_correlation_time_rates);

    EnzymeConcentrations.remove(0)

    EnzymeConcentrations = pd.DataFrame(EnzymeConcentrations);

    EnzymeConcentrations = EnzymeConcentrations.T;

    EnzymeConcentrations.columns = high_correlation_time_rates.columns;

    high_corrs_and_EZCs = [EnzymeConcentrations, high_correlation_time_rates]

    normalized_high_correlation_time_rates = pd.concat(high_corrs_and_EZCs, axis = 0, ignore_index = False);

    normals = normalized_high_correlation_time_rates.iloc[1:, :] / normalized_high_correlation_time_rates.iloc[0, :];
    
    linearized_rates_2 = normals.mean(axis = 1);    

    linear_rate_stdevs = normals.std(axis = 1);

    high_correlation_time = pd.DataFrame(high_correlation_time);

    high_correlation_time.index = linearized_rates_2.index

    max_lin_plot_data = [high_correlation_time, linearized_rates_2]

    max_lin_plot_data = pd.concat(max_lin_plot_data, axis = 1, ignore_index = True);

    max_linearized_rate_plot = max_lin_plot_data.plot.bar(x = 0, yerr = linear_rate_stdevs, legend = None);

    max_linearized_rate_plot.set_xlabel('Timepoints with high linearity between protein loading & reaction rate');

    max_linearized_rate_plot.set_ylabel('Rate of reaction mM/min');

    max_linearized_rate_plot.set_title('Maximum rates at linearized timepoints');
    
    return time, mean_data, rate_times, rate_zeros, linearity, rsq_plot_vals, max_lin_plot_data, linear_rate_stdevs 












