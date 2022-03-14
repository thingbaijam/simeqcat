import random 
from datetime import datetime 
import time


def get_random_dates(nevents,\
                        start_date, end_date):
    event_dates = []
    for i in range(nevents):
        while True:
            r = random.random()
            if (r>0.0) & (r<1.0):
                break;
                
        tdate = start_date + (end_date - start_date) * r
        event_dates.append(tdate)
    event_dates.sort()
    return event_dates
 

def random_datetime(min_year=1900, max_year=datetime.now().year):
    # generate a datetime in format yyyy-mm-dd hh:mm:ss.000000
    start_year = datetime(min_year, 1, 1, 0, 0, 0)
    end_year = datetime(max_year, 12, 31, 23, 59, 0)
    return start_year + (end_year - start_year) * random.random()




def generate_random_datetimes(min_year=1900, \
                              max_year = datetime.now().year,\
                              nevents=1):
    
    event_dates = []
    for i in range(nevents):
        event_dates.append(random_datetime(min_year=min_year, max_year=max_year))

    event_dates.sort()
    return event_dates
    
    

def get_yearfraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction
