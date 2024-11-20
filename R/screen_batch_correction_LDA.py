# -*- coding: utf-8 -*-
"""

@author: arshia hassan
"""

import argparse
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

'''
Function:
    filter_screen()
Arguments:
    screen_list : List of compound screen-labels (label format example: CHEM014_Bortezomib_T18 ) (required)
    samples_in_batch: Minimum number of screens to form a batch (default = 2)
Return: 
    List of labels that meet the minimum criteria.
Description:
    Objective: Extract/filter the screen-labels that meet the minimum batch size criteria
    - Needed in the Batch correction step
    
'''
def filter_screen(screen_list,samples_in_batch = 2):
    #create dataframe to store label (CHEM014_Bortezomib_T18), screen-batch number (CHEM014), compound name (Bortezomib) and time (T18)
    screen_info = pd.DataFrame(columns=['labels','screen', 'chemical','time'])
    rowIndex = 0
    #for each screen label extract the sub-parts and store in screen dataframe
    for label in screen_list:        
        temp_screen,temp_chemical,temp_time = label.split('_')  
        screen_info.loc[rowIndex, 'labels'] = str(label)
        screen_info.loc[rowIndex, 'screen'] = str(temp_screen)
        screen_info.loc[rowIndex, 'chemical'] = str(temp_chemical.upper())
        screen_info.loc[rowIndex, 'time'] = str(temp_time)
        rowIndex+=1
    # Count the number of screens in each screen-batch based on screen-batch number (i.e. CHEM014) 
    screen_count = pd.DataFrame(columns=['screen', 'count'])
    screen_count['screen'] = screen_info.screen.value_counts().index.tolist()
    screen_count['count'] = screen_info.screen.value_counts().tolist()
    
    #Extract the screen-batch numbers that meet the minimum batch size criteria
    screen_greater = list(screen_count[screen_count['count'] > samples_in_batch]['screen'])
    screen_labels = list(screen_info[screen_info['screen'].isin(screen_greater)]['labels'])

    return screen_labels

'''
Function:
    create_binary_standard()
Arguments:
    screen_labels : (required)
    batch_dictionary : Batch-screen and batch-label mapping (required)
    set_na (type:binary): If True, same compunds in the same batch are set to NA - don't penalize them for being in the same screen-batch, because they should be in the same batch(default = True)
Return: 
    1D array with binary standard 
Description:
    Objective: Create binary standard for screen-label pairs 
    (1:if they belong to the same screen-batch , 0:if they do not belong to the same screen-batch )
    - Needed to generate ROC score
    
'''
def create_binary_standard(screen_labels, batch_dictionary, set_na = True):
    screens_length = len(screen_labels) # number of screen labels
    binary_matrix = np.zeros((screens_length, screens_length)) # Initiate pair-wise binary label matrix to 0
    for screen_index_1 in range(0,screens_length): # for 1st screen in screen-label-pair
        screen_1 = screen_labels[screen_index_1] # get screen label
        screen_1_list = screen_1.split('_') # extract screen information 
        for screen_index_2 in range(0,screens_length): # for 2nd screen in screen-label-pair
            screen_2 = screen_labels[screen_index_2] # get screen label
            screen_2_list = screen_2.split('_') # extract screen information 
            if(batch_dictionary[screen_1]==batch_dictionary[screen_2]) : # if both screens are from the same batch as depicted in batch_dictionary               
                if set_na and (screen_1_list[1] == screen_2_list[1]): # if both screens are from the same compound
                    binary_matrix[screen_index_1][screen_index_2] = np.nan # set to NA (We don't want to penalize same compunds for batch effect)
                else: # otherwise
                    binary_matrix[screen_index_1][screen_index_2] = 1 # set to 1 (the screen-label-pair is in the same screen batch)
    binary_1d = binary_matrix.flatten() # Transform matrix to 1D array
    binary_1d = binary_1d[~np.isnan(binary_1d)] # remove NA entries    
    return binary_1d

'''
Function:
    create_pcc_scores()
Arguments:
    screen_data : gene profiles across compounds (gene x screen)(required)
    batch_dictionary : Batch-screen and batch-label mapping (required)
    set_na (type:binary): If True, same compunds in the same batch are set to NA - don't penalize them for being in the same screen-batch, because they should be in the same batch(default = True)
Return: 
    1D array with Pearson Correlation Co-efficient scores
Description:
    Objective: Generate PCC scores for screen-label pairs 
    - Needed to generate ROC score
    
'''
def create_pcc_scores(screen_data, batch_dictionary, set_na = True):
    screen_labels = screen_data.columns.values # Get screen labels
    screens_length = len(screen_labels) # number of screen labels 
    data_t = screen_data.T # transpose data (screen x genes)
    corr_matrix = np.corrcoef(data_t)  # calculate Pearson Correlation Coefficient among screen profiles 
    if set_na:
        for screen_index_1 in range(0,screens_length): # for 1st screen in screen-label-pair
            screen_1 = screen_labels[screen_index_1] # get screen label
            screen_1_list = screen_1.split('_') # extract screen information
            for screen_index_2 in range(0,screens_length): # for 2nd screen in screen-label-pair
                screen_2 = screen_labels[screen_index_2] # get screen label
                screen_2_list = screen_2.split('_') # extract screen information
                if(batch_dictionary[screen_1]==batch_dictionary[screen_2]) : # if both screens are from the same batch as depicted in batch_dictionary
                    if(screen_1_list[1] == screen_2_list[1]): # if both screens are from the same compound
                        corr_matrix[screen_index_1][screen_index_2] = np.nan  # set to NA (We don't want to penalize same compunds for batch effect)                  
    corr_1d = corr_matrix.flatten() # Transform matrix to 1D array
    corr_1d =corr_1d[~np.isnan(corr_1d)]  # remove NA entries 
    return corr_1d

'''
Function:
    rpc_plot_generation()
Arguments:
    fpr: false positive rate calculated by roc_curve() function
    tpr: true positive rate calculated by  roc_curve() function
    roc_auc_score_: rocauc score calculated by roc_auc_score()
    file_path : output plot file path
Description:
    Generate ROC plot
'''
def rpc_plot_generation(fpr,tpr,roc_auc_score_,file_path):
    plt.figure()
    lw = 2
    plt.plot(
        fpr,
        tpr,
        color="darkorange",
        lw=lw,
        label="ROC curve (area = %0.5f)" % roc_auc_score_,
    )
    plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("")
    plt.legend(loc="lower right")
    plt.savefig(file_path)
    plt.close()
    return

'''
Function:
    run_batch_correction()
Arguments:
    input_file_path: Path to file to apply LDA batch correction
    output_file_directory: Directory to save roc plot files
Return:
    Path of the output file
Description:
    Applies LDA batch correction to data and removes LD components until the ROCAUC score drops below .51
   
'''
def run_batch_correction(data, output_file_directory):
    print(str('Running batch correction'))
    
    if "gene" in data.index.values: 
        data.set_index("gene", inplace=True) # set gene column as index
        
    screen_list = data.columns.values # get all screen labels
    
    #extract screen and replicate batches
    screen_labels_batch = filter_screen(screen_list, samples_in_batch = 1) # Extract batches with at least two screens  #updated
    screen_orphan_batch = list(set(data.columns).difference(screen_labels_batch)) # Rest of the screens make up Orphan batch
    
    #set up screen-label to batch mapping for screen batches
    screen_batch_dictionary = {}
    orphan_dictionary = { }
    for label in screen_list:
        temp_screen,temp_chem,temp_time = label.split('_')
        if label in screen_labels_batch: #
            screen_batch_dictionary[label] =  temp_screen
        elif label in screen_orphan_batch:
            orphan_dictionary[label] =  "CHEMORPHAN"
    screen_batch_dictionary.update(orphan_dictionary) 
    
    #create dataframe with batches
    multiple_batch_df = data[screen_labels_batch]
    orphan_batch_df = data[screen_orphan_batch]
    frames = [multiple_batch_df,orphan_batch_df]
    df = pd.concat(frames,axis = 1)
    
    print(str('Training LDA'))
    # LDA analysis
    X = df.T # Shape: (batch, genes)
    y = np.array(list(screen_batch_dictionary.values()))
    lda_eigen = LinearDiscriminantAnalysis(solver='eigen', n_components=9, shrinkage='auto', store_covariance=True)
    lda_eigen.fit(X,y)
    #X_transform_eigen = lda_eigen.transform(X)
    print(str('Training LDA complete.'))
    
    # scaling matrix: Scaling of the features in the space spanned by the class centroids. Shape: (gene, gene)
    pc_matrix_temp = lda_eigen.scalings_
    pc_matrix_dup = np.copy(pc_matrix_temp)
    row_l2_norm_L =  np.linalg.norm(pc_matrix_dup, axis=1)
    pc_matrix_dup /= row_l2_norm_L[:, None]
    
    main_matrix = df.T.to_numpy() # Shape: (batch, genes)
    main_matrix= main_matrix.astype('float')
    
    components_reduced = 0
    #Remove components
    while components_reduced > -1: 
        print(str(components_reduced))
        
        # L shape: (gene,component)
        L = pc_matrix_dup[:,0:components_reduced] # components_reduced
        # Lt shape: (component,gene)
        Lt = pc_matrix_dup[:,0:components_reduced].T 
    
        # Reduced matrix shape (batches, genes)
        reduced_matrix = np.linalg.multi_dot([main_matrix,L,Lt])
    
        # Removed PC matrix: shape (batches, genes)
        pc_removed_matrix = main_matrix - reduced_matrix    
        pc_removed_matrix_df = pd.DataFrame(pc_removed_matrix)
        pc_removed_matrix_df.columns = X.columns.values
        pc_removed_matrix_df.index = X.index.values
        pc_removed_matrix_df = pc_removed_matrix_df.T
        
        # calculate roc_auc score to evaluate batch remove effect, we want a score close to random
        # screens in the same batch shouldn't be more correlated than random screen pairs
        pc_removed_screen_list = pc_removed_matrix_df.columns.values        
        screen_binary_na = create_binary_standard(pc_removed_screen_list, screen_batch_dictionary, set_na = True) # create binary scores (ture)
        screen_pcc_na = create_pcc_scores(pc_removed_matrix_df, screen_batch_dictionary, set_na = True) # create PCC scores (predicted)
        roc_auc_score_screen_pcc_na = roc_auc_score(screen_binary_na, screen_pcc_na) # calculate ROCAUC
        print(str(components_reduced) + '\tscreen auc dc\t' + str(roc_auc_score_screen_pcc_na))
        
        filepath = output_file_directory + 'bc_lda_'+ str(components_reduced)
        fpr, tpr, _ = roc_curve(screen_binary_na, screen_pcc_na)
        rpc_plot_generation(fpr,tpr,roc_auc_score_screen_pcc_na,filepath+'.png')
        
        #if roc_auc score drops below .51 save data file and plot, stop batch correction        
        if np.around(roc_auc_score_screen_pcc_na,decimals=2) < 0.51 or components_reduced > 20: 
            print('LDA correction complete.')
            print('Removed components: ' + str(components_reduced))
            
            break
        components_reduced = components_reduced + 1
    return pc_removed_matrix_df
