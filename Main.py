import lib_GLONASS as brd

year = 'XXXX'   #Insert the year of interest 
day = 'DD' #Doy of year
day = str(i).zfill(3)
dir = 'Destination of your data'

file = dir + 'brdm'+day+'0.'+str(year)[2:4]+'p' # Opening the broadcast file 
fl = open(file, 'r' )

cor = open(dir + '/resGLONASS'+year+'_'+day+'.csv', 'w') #creating file for residuals
fl = fl.readlines()

#Calling the function that will calculate the residuals
brd.CalcCoordGLONASS(fl,cor,dir)


