#!/usr/bin/env python
"""
Generate a random sequence with a seed so we know what the sequence was to
use to spinup GDAY

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (20.06.2015)"
__email__ = "mdekauwe@gmail.com"

import numpy as np
import calendar
import sys
import datetime as dt 

def main():
    
    # Set the seed so we can repeat this if required
    np.random.seed(42)
    
    start_yr = 1960 # weird stuff in the 1950s AWAP data
    end_yr = 1990
    out_yrs = end_yr - start_yr
    yrs = np.arange(start_yr, end_yr+1)
    
    # preserve leap yrs, so find them first
    leapyrs = np.zeros(0)
    for yr in yrs:
        if calendar.isleap(yr):    
            leapyrs = np.append(leapyrs, yr)
        #print >>f, float(yr) 
    
    # However we don't want the leapyrs in the sequence, so exclude them
    yrs = np.array([yrs[i] for i, yr in enumerate(yrs) if yr not in leapyrs])
    
    shuff_years = randomise_array(out_yrs, yrs)
    shuff_years_leap = randomise_array(out_yrs, leapyrs)
    
    
    i = 0
    for yr in np.arange(start_yr, end_yr):
        
        if i == 0:
            prev_yr_leap = 1666 # anything not in the sequence
            prev_yr = 1666 # anything not in the sequence
            
        if calendar.isleap(yr):
            out_yr = shuff_years_leap[i]
           
            # Make sure we don't get the same year twice
            while prev_yr_leap == int(out_yr):
                i += 1
                out_yr = shuff_years_leap[i]
                
            print "%d," % (int(out_yr)),
            prev_yr_leap = shuff_years_leap[i]
        else:
            out_yr = shuff_years[i]
            
            # Make sure we don't get the same year twice
            while prev_yr == int(out_yr):
                i += 1
                out_yr = shuff_years[i]
            
            print "%d," % (int(out_yr)),
            prev_yr = shuff_years[i]
        
        i += 1
    print
       
def randomise_array(out_yrs, yrs):
    # make a sequence longer than the number of years we actually want
    num_seq = np.ones(out_yrs * len(yrs))
    num_years = len(yrs) * len(num_seq)
    shuff_years = (yrs * num_seq[:,None]).reshape(num_years)
    np.random.shuffle(shuff_years)
    
    return shuff_years
        
   
if __name__ == "__main__":

    main()
    