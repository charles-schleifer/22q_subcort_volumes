#!/usr/bin/env python
# coding: utf-8

# import abagen and other libraries
import os
# pip install abagen
import abagen
from abagen import images
from abagen import reporting
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)

# Download AHBA data
# TODO: uncomment if you need to re-download the data
#files = abagen.fetch_microarray(donors='all', verbose=1)

# get path to parent directory (github repo)
# this works if full script is run, otherwise, need to set project manually to github repo path
# project='/Users/charlie/Dropbox/github/22q_spatial_transcriptomic/'
project=os.path.dirname(__file__)+'/'    

# example data
exatlas = abagen.fetch_desikan_killiany(native=True)
#/Users/charlie/opt/anaconda3/lib/python3.9/site-packages/abagen/data/

# get freesurfer data for all donors
#/Users/charlie/abagen-data/freesurfer
fs = abagen.datasets.fetch_freesurfer(donors=('9861', '10021', '12876', '14380', '15496', '15697'))

# create new atlas dictionary for native space individual thalamus
fspath='/Users/charlie/abagen-data/freesurfer/' 
thal=exatlas
# set paths in dictionary to thalamus image
for i in thal['image']:
    thalpath=(fspath+'/'+'donor'+i+'/mri/ThalamicNuclei.nii.gz')
    thal['image'][i]=thalpath
# manually make info file for thalamus from freesurfer LUT and add to atlas
thal['info']=project+'thal_info_abagen.csv'


# volume only
#ahbaThalNative = abagen.get_expression_data(atlas=thal['image'], atlas_info=thal['info'], norm_structures=True, return_donors=True,donor_probes='independent', verbose=2)
ahbaThalNative = abagen.get_expression_data(atlas=thal['image'], atlas_info=thal['info'], norm_structures=True, return_donors=True,lr_mirror='bidirectional', verbose=2)


# save results
#for d in ahbaThalNative.keys():
#    ahbaThalNative[d].to_csv(project+'/thal_native_expression.'+d+'.csv')
    

# generate and save methods report
#volReport = reporting.Report(atlas=thal['image'], atlas_info=thal['info'], norm_structures=True, return_donors=True,donor_probes='independent')
#volReport_out = volReport.gen_report()
#with open(project+'/abagen_methods_report_thal_native.txt', 'w') as text_file:
#    text_file.write(volReport_out)


# get regions where all donors have data
regions={}
for d in ahbaThalNative.keys():
    regions[d]=set(ahbaThalNative[d].index[ahbaThalNative[d].notna().any(axis=1)])
    
    
regions_in_all=regions['9861'].intersection(regions['10021'],regions['12876'],regions['14380'],regions['15496'],regions['15697'])

print(regions_in_all)






    