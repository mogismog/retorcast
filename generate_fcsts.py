#!/usr/bin/env python
#
# Python program for outputting probabilistic tornado forecasts using an analog method with
# the 2nd Generation Reforecast Dataset (see: Hamill, et al. 2013 BAMS article).
# 
# You can run the program just inputting the
# date in YYYYMMDD00 format, a lead time (in days) and an ROI (80/160/240). For example, if we want
# a Day 4 forecast from January 21, 2013, we'd write:
#
#   python retorcast.py 2013012100 4
#
# ****************************IMPORTANT!************************************
# Be sure to change the configuration file (config.ini) to appropriate settings!
#
# Requires numpy, py-netCDF4, pygrib, retorcast modules(should be included...),
# and assorted (in folder ./data/).
#
# Version: 1.0 (1/2014) - Initial version, uses analog method, super clean and easy!
#
# Contact:  Francisco M. Alvarez <falvare6@slu.edu>
#
# copyright (c) by Francisco M. Alvarez.
# 
# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose and without fee is hereby granted,
# provided that the above copyright notice appear in all copies and that
# both the copyright notice and this permission notice appear in
# supporting documentation.
# THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
# EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
#

import sys
from datetime import datetime
from retorcast.analog import generate_probabilities,multi_plot_analog

if len(sys.argv) < 3 or len(sys.argv) > 3:
    raise Exception('\nUsage:\n python generate_fcsts.py YYYYMMDD00 lead_time_in_days\n\nTo run all forecasts (Days 1-10) for the current date, use:\n python generate_fcsts.py today all')
    
# --- parsing command line arguents!
try:
    lead_time = int(sys.argv[2]) # --- lead time (days)
except ValueError:
    lead_time = sys.argv[2]
    
fdate = sys.argv[1] # --- Forecast date
    
# --- Turning forecast date into datetime object:
if fdate.lower() == "today":
    fdate = datetime.utcnow()
    fdate = fdate.replace(hour=0,minute=0,second=0,microsecond=0)
else:
    fdate = datetime(int(fdate[:4]),int(fdate[4:6]),int(fdate[6:8]))

# --- Call the main function to generate analogs...
try:
    if lead_time.lower() == "all":
        for lead_date in xrange(1,11,1):
            generate_probabilities(fdate,lead_date)
            #multi_plot_analog(fdate,lead_date)
            
except AttributeError:
    #generate_probabilities(fdate,lead_time,use_pct=True)
    multi_plot_analog(fdate,lead_time)
