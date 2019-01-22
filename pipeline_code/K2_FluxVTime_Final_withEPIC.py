# -*- coding: utf-8 -*-
"""
@author: Sheila Sagear
Used EPIC ID
Flux scaling? - normalized to 1
***ALLWISE/EPIC/K2 ID? Some coordinates have multiple targets, some have none***
Time and Flux units? (relative brightness)
Is flux comparable to brightness for these purposes?
"""
#ymax as y limit
from pylab import xlabel, ylabel, title, plot, show, ylim

#count number of files being graphed
filecount = sum(1 for line in open('all_txt_files.txt'))

#opens list of file names, creates list of file names
with open('all_txt_files.txt') as f:
    content = f.readlines()

#converts each file name to string, removes \n from end, opens each object's txt file (loop)
#and stores (x,y) in list 'data'
tempindex = 0
while (tempindex < filecount):
    file = str(content[tempindex])
    file = file.strip('\n')
    file = file.strip('\r')
    with open(file) as f:
        data1 = f.read()

#converts txt file to string (and removes heading: 30 characters) 
#in order to isolate each (x,y) coordinate as an element in list 'data'
        datastr = str(data1)
        datastr = datastr[30:]
        data1 = datastr.split('\n')


#removes comma after each (x,y) coordinate;
#isolates x and y values as indicies of list 'data'
        index = -1
        while (index < len(data1)):
            tempstring = str(data1[index])
            data1[index] = tempstring.rstrip(',')
            data1[index] = tempstring.split(', ')
            index+=1
        data1.pop
        
        data2 = sum(data1, [])
        
        index = 0
        while (index < len(data2)):
            if index % 2 == 1:
                data2[index] = data2[index].rstrip(',')
            index+=1
        data2.pop()
           
#converts str data points to float
        data_final = [float(i) for i in data2]
    
#defines x and y values by index of 'data'
        x_list = data_final[0::2]
        y_list = data_final[1::2]

#normalizes flux values to 1.0 (at avg of flux values)
        y_count = len(y_list)
        sumy_list = sum(y_list)
        y_avg = sumy_list/y_count
        
        y_list = [i/y_avg for i in y_list]
        
        xlabel('Time')
        ylabel('Corrected Flux (normalized to 1.0)')

#targets titled by EPIC names
        targetname = file[25:]
        targetname = targetname[:-35]
        title('EPIC ' + targetname + ' Light Curve')
#not sure if y limits are appropriate
        #ylim([.9, 1.5])
        plot(x_list,y_list)
        ylim(-.005, 2.005)
        show()

    
    tempindex+=1