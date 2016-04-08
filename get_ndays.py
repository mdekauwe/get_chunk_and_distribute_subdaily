import calendar


year_seq = [1964, 1962, 1970, 1989, 1968, 1985, 1973, 1977, 1972,\
            1969, 1987, 1983, 1984, 1982, 1962, 1971, 1968, 1962,\
            1965, 1990, 1960, 1977, 1969, 1966, 1968, 1965, 1982,\
            1985, 1980, 1966]
ndays = 0
for yr in year_seq:
    if calendar.isleap(yr):
        ndays += 366
    else:
        ndays += 365

print ndays

ndays = 0
for yr in xrange(1990, 2011+1):
    if calendar.isleap(yr):
        ndays += 366
    else:
        ndays += 365

print ndays
