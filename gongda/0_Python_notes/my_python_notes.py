#####################################################################################################
# my incomplete notes on Python basics

# string, list, dictionary

# Flow control
# for loops - list comprehension 
# if else statements
# while - continue - break
# range()
# with statements
#####################################################################################################
# string

# 1> access characters or slice the string

str = "Leicester Supersite NOx"

print(str)
print(str[0],str[-1],sep = ",")
print(str[0:9])
print(str[:-3])
print(str[0:9],str[20:23])
print(type(str),len(str),end=' ') # print results in one line 
print("L" in str) # check if element "L" is in the string               

# 2> string.split(separator, maxsplit) 

DateTime = "2016-01-01 12:59:59"
Date  = DateTime.split(" ")[0]
Time  = DateTime.split(" ")[1]
Year  = Date.split("-")[0]
Month = Date.split("-")[1]
Day   = Date.split("-")[2]
Hour  = Time.split(":")[0]
Min   = Time.split(":")[1]
Sec   = Time.split(":")[2]

print(Year,Month,Day,Hour,Min,Sec, sep = "   ")

# reverse the process above
date = [Year,Month,Day]
sep = '-'
date = sep.join(date)

time = [Hour,Min,Sec]
sep = ':'
time = sep.join(time)

datetime = date + " " + time
print(datetime)

# 3> replace the string
str_1 = "Leicester Supersite NOx"
str_2 = str_1.replace("NOx","NO2")
print(str_2)
#####################################################################################################
# list 

# example 0
surface_NO2 = [] # create an empty list 

# example 1
a = [1,2,3,4]
a[2] = "x"       # assign a[2]
b = a + [5,6,7]
c = [x*2 for x in b]

for x in [a,b,c]:
    print(x)
    
# example 2
a = ['Leicester','Cambridge','London','York','Edinburgh']
print(a[2:4])
print(a[-2:])
a.insert(4,'Nottingham') # add "Nottingham" as a[4]
del a[0:2]               # delete a[0], a[1]
a.remove('Edinburgh')    # remove a certain element
a.append('Leeds')        # add it automatically to the end
print(a)

# example 3
a = [1,2,5,6,3,2,2,0]
b = sorted(a)         # sorted() creates a new list while keeping the original one
b.count(2)            # count the occurrences of element "2"
a.sort()              # .sort() changes the original list        
#####################################################################################################
# dictionary 

# structure: {key: value} "key" must be Integer or String - "value" can be anything

# example 1
gas = {'NO2': 'nitrogen dioxide',
       'SO2': 'sulfur dioxide'}

print(gas)
print(gas.keys())
print(gas.values())

# example 2
short_name = ['NO2','SO2','CO','O3','PM2.5']
long_name  = ['nitrogen dioxide', 'sulfur dioxide', 'carbon monoxide', 'ozone', 'particulate matter 2.5']

pollutants = {}

for i in range(len(short_name)):
    dict1 = "{}".format(short_name[i])
    dict2 = "{}".format(long_name[i])
    pollutants[dict1]=dict2
    
print(pollutants)

pollutants['HCHO'] = 'formaldehyde' # add one more pollutant
del pollutants['NO2']               # remove NO2
print(pollutants['SO2'])            # look up SO2 
print('SO2' in pollutants)          # check if SO2 is in dictionary
print(pollutants)                   
#####################################################################################################
#####################################################################################################
# Flow control
# for loops - list comprehension
# if else
# while - continue - break

# example 1
a = [1,2,3,4,5,6,7,8]
b = []
c = []

for i in range(len(a)):
    if a[i] > 3 and a[i] < 6:
        b.append(a[i])
    else:
        b.append('NA')

for i in range(len(a)):
    if a[i] > 3 or a[i] < 6:
        c.append(a[i])
    else:
        c.append('NA')

print(b)
print(c)

# if all([condition_1,..., condition_n):
# if any([condition_1,..., condition_n):

# example 2
a = [-2,-1,0,1,2]

b = []
for i in range(len(a)):
    if a[i] < -1:
        b.append('too low')
    elif -1 <= a[i] < 0:
        b.append('low')
    elif a[i] == 0:
        b.append('good')
    elif 1 <= a[i] < 2:
        b.append('high')
    else:
        b.append('too high')

c = ['too low' if i < -1 else 'low' if -1 <= i < 0 
                         else 'good' if i == 0 
                         else 'high' if 1 <= i <2
                         else 'too high'for i in a]

d = [i*100 if i < -1 else i*200 if -1 <= i < 0 
                     else i*300 if i == 0 
                     else i*200 if 1 <= i <2
                     else i*100 for i in a]

for x in [a,b,c,d]:
    print(x)

# example 3

i = 0
while i < 5:
    print(i, end = "  ")
    i = i+1

for i in range(15):
    if i%5 == 0:
        continue
    print(i, end = "/ ")

for i in range(15):
    print(i, end = "  ")
    if i >= 5:
        break        
#####################################################################################################
# other notes that might be useful

==   # equal to
!=   # not equal to
<=   # less than or equal to
>=   # more than or equal to

A and B # "&"
A or B  # "|"
not A  

a in A
a not in A

range(n)      # 0,1,2...n-1
range(1,n)    # 1,2,3...n-1

# with statement
with xxx as f:
    data = f.xxx
#####################################################################################################
# other tricks that might be useful

# Trick 1: format a string
data = 'TROPOMI NO2'
date = '01/01/2020'

print("{0} is great".format(data))
print("{0} ({1})".format(data,date))
print("{0}\n{1}".format(data,date))

# Trick 2: print results seperately to make it easier to read in Jupyternotebook
test = Dataset(TropOMI_files[0], "r", format="NETCDF4") 
print(test,test.groups['PRODUCT'],test.groups['METADATA'], sep = '\n##########################\n')

# End
#####################################################################################################
