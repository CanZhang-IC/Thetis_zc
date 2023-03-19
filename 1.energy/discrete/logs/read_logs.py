import  numpy as np
# for name in ['rec','sta','all','noyaw-rec','noyaw-sta']:

f = open('noyaw-sta.log','r')
lines = f.readlines()
power_all_time = []
for line in lines:
    if 'Current and integrated power for each turbine' in line :
        a = line.index('[')
        b = line.index(']')
        power_one_time = line[a+2:b-1].split(',')
        power_one_time_float = []
        for i in power_one_time:
            power_one_time_float.append(float(i))
        power_all_time.append(power_one_time_float)

# power_array = np.array(power_all_time)
# power_average1 = np.average(power_array,axis=0)
# power_average2 = np.average(power_array,axis=1)
# print(power_average1,power_average2)

aa =[15059938.295935603, 6461910.710982086, 5070302.496777756, 40438305.29598973, 36123182.63422021, 20737026.222772222, 23795045.48594738, 31589508.424124565, 41136421.64255234, 10499257.915163139, 12761234.806500144, 16863982.352332402, 5368372.8694443805, 5072600.204186277, 6376516.227242672, 4782799.020218644, 2399865.111986514, 4235435.381217366]
aa = np.array(aa)/46800
print(sum(aa))