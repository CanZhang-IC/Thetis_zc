
for i in range(20):
    i += 1
    filename = './5min-1core-sediment-exner'+str(float(i))+'/log'
    with open(filename,'r',encoding='utf-8') as f :
        lines = f.readlines()
        the_line = lines[-3]
        print(i,the_line[:5])