from datetime import datetime,timedelta
import itertools

def get_1mo_dates(inyr,inmo,indate,byear,eyear):
    """
 Used with netCDF4 files and py-netCDF.

 From 1985-2011, a number of days range to search
 for dates centered around the given date and forecast time,
 returns the applicable dates in a list of daetime objects.
 In other words, returns a list of file name dates for us to
 search for analogs/use in logistic regression/whatever.

 indate - Initial date (1,31)
 inyr - Initial year, YYYY (1985 - )
 inmo - Initial month, (1,12)
 window - range of dates in past years to search, e.g. 45 will find dates 45 days before/after indate

 Returns:
 outdates - List of dates meeting the criteria
    """

    fnlist = []

    #print inmo,indate
    try:
        xdate = datetime(byear,inmo,indate)
    except ValueError:
        xdate = datetime(byear,inmo,indate-1)
    else:
        xdate = datetime(byear,inmo,indate)
    while xdate < datetime(eyear+1,1,1):
        #print xdate
        if xdate.year == inyr:
            try:
                xdate = datetime((xdate.year + 1),inmo,indate)
            except ValueError:
                xdate = datetime((xdate.year + 1),inmo,indate-1)
            continue
        for datechange in xrange(0,35):
                tdelta = timedelta(days=datechange)
                analogdate = xdate + tdelta
                #print analogdate,xdate
                if analogdate.year > eyear:
                    continue
                if analogdate.year == inyr:
                    continue
                if analogdate.month != xdate.month:
                    continue
                fnlist.append(analogdate)
        try:
            xdate = datetime((xdate.year + 1),inmo,indate)
        except ValueError:
            xdate = datetime(xdate.year+1,inmo,indate-1)

    return fnlist

def get_analog_dates(inyr,inmo,indate,byear,eyear,**kwargs):
    """
 Used with netCDF4 files and py-netCDF.

 From 1985-2011, a number of days range to search
 for dates centered around the given date and forecast time,
 returns the applicable dates in a list of daetime objects.
 In other words, returns a list of file name dates for us to
 search for analogs/use in logistic regression/whatever.

 indate - Initial date (1,31)
 inyr - Initial year, YYYY (1985 - )
 inmo - Initial month, (1,12)
 window - range of dates in past years to search, e.g. 45 will find dates 45 days before/after indate
 byear - earliest year for potential dates (usually 1985)
 eyear - latest year for potential dates (usually 2011)

 Returns:
 outdates - List of dates meeting the criteria
    """

    bias_corr = kwargs.get('bias_corr',False)


    fnlist = []
    date_list = []
    #print inmo,indate
    try:
        xdate = datetime(byear,inmo,indate)
    except ValueError:
        xdate = datetime(byear,inmo,indate-1)
    else:
        xdate = datetime(byear,inmo,indate)
    while xdate < datetime(eyear+1,1,1):
        #print xdate
        if xdate.year == inyr:
            try:
                xdate = datetime((xdate.year + 1),inmo,indate)
            except ValueError:
                xdate = datetime((xdate.year + 1),inmo,indate-1)
            continue
        for datechange in reversed(xrange(0,100)):
            if xdate.month > 1:
                tdelta = timedelta(days=datechange)
                analogdate = xdate - tdelta
                if analogdate.year < byear:
                    continue
                if analogdate.year == inyr:
                    continue
                if analogdate < datetime(xdate.year,xdate.month-1,1):
                    continue
                fnlist.append(analogdate)
            elif xdate.month == 1:
                tdelta = timedelta(days=datechange)
                analogdate = xdate - tdelta
                if analogdate.year < byear:
                    continue
                if analogdate.year == inyr:
                    continue
                if analogdate < datetime(xdate.year-1,12,1):
                    continue
                fnlist.append(analogdate)
        for datechange in xrange(1,101):
            if xdate.month < 12:
                tdelta = timedelta(days=datechange)
                analogdate = xdate + tdelta
                if analogdate.year > eyear:
                    continue
                if analogdate.year == inyr:
                    continue
                try:
                    datetime(xdate.year,xdate.month+2,1)
                except ValueError: # --- xdate.month == 11
                    if analogdate >= datetime(xdate.year+1,1,1):
                        continue
                else:
                    if analogdate >= datetime(xdate.year,xdate.month+2,1):
                        continue
                fnlist.append(analogdate)
            elif xdate.month == 12:
                tdelta = timedelta(days=datechange)
                analogdate = xdate + tdelta
                if analogdate.year > eyear:
                    continue
                if analogdate.year == inyr:
                    continue
                if analogdate >= datetime(xdate.year+1,2,1):
                    continue
                fnlist.append(analogdate)
        try:
            xdate = datetime((xdate.year + 1),inmo,indate)
        except ValueError:
            xdate = datetime(xdate.year+1,inmo,indate-1)


    # --- Here, since we are now using bias-corrected data, we can get additional potetial analog dates!
    if bias_corr:

        date_list.append(fnlist)

        #for n_mo in xrange(1,13,1):
        #    if (n_mo >= inmo-1) and (n_mo <= inmo+1):
        #        continue
        #    else:
        #        date_list.append(get_1mo_dates(int(inyr),n_mo,1,byear,eyear))
        if (inmo < 2) or (inmo > 9):
           date_list.append(get_1mo_dates(int(inyr),3,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),4,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),5,1,byear,eyear))
        if (inmo == 2):
           date_list.append(get_1mo_dates(int(inyr),4,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),5,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),10,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),11,1,byear,eyear))
        if (inmo == 3):
           date_list.append(get_1mo_dates(int(inyr),5,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),10,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),11,1,byear,eyear))
        if (inmo == 4):
           date_list.append(get_1mo_dates(int(inyr),9,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),10,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),11,1,byear,eyear))

        # --- Now flatten and return the list
        date_list = list(itertools.chain.from_iterable(date_list))
        return date_list
    else:
        return fnlist