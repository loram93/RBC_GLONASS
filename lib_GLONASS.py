import math
import datetime
from numpy.linalg import inv
from numpy import *

def CalcCoordGLONASS(fl, cor,dir):
    for lines in fl:
        if 'LEAP SECONDS' in lines:
            lps = int(lines[0:10])
    s,hd = sats2(fl)
    for lines in fl:
        if lines.startswith('G01'):
            line = lines.split()
            day = int(line[3])
            month = int(line[2])
            yr = int(line[1])
            year = yr
            break
    doy = datetime.datetime(year, month,day).timetuple().tm_yday
    dow = datetime.datetime(year, month,day).weekday() +1
    if dow == 7:
        dow = 0
    #print(dow)
    d = str(doy).zfill(3)
    gl = open(dir + str(year)+'_'+d, 'r')
    gl = gl.readlines()



    print(day, month, year)
    for sat in s:
        for j in range(hd,len(fl)):
            sc = 0
            if fl[j].startswith(sat):
                if sat[0] in ['R']:
                    # print sat
                    # Constants #########################################
                    mu = 398600.44
                    c20 = -1082.6257e-6
                    ae = 6378.136
                    we = 0.7292115e-4
                    Tu = jdn(yr, month, day) / 36525
                    H0 = (24110.54841 + 8640184.812866 * Tu + 0.093104 * (Tu ** 2) - 6.2e-6 * (
                            Tu ** 3))  # *(2*pi)/86400
                    w3 = (we + 4.3e-15 * Tu) * 86400 / (2 * pi)

                    ####################################################
                    ## Elements coming from the the nav file          ##
                    ####################################################
                    TauN = float(fl[j][23:42])
                    GammaN = float(fl[j][42:61])
                    tb = int(float(fl[j][61:80]))
                    X0 = float(fl[j + 1][0:23])
                    Xdot = float(fl[j + 1][23:42])
                    xdotdot = float(fl[j + 1][42:61])
                    health = float(fl[j + 1][61:70])
                    Y0 = float(fl[j + 2][0:23])
                    Ydot = float(fl[j + 2][23:42])
                    ydotdot = float(fl[j + 2][42:61])
                    Z0 = float(fl[j + 3][0:23])
                    Zdot = float(fl[j + 3][23:42])
                    zdotdot = float(fl[j + 3][42:61])


                    if X0 == 0:
                        break

                    ## Start of calculations #############################
                    ######################################################

                    for i in range(tb - 15 * 60, tb + 15 * 60):
                        if (i % 900 == 0):
                            Ste = (H0 + (we * (i - 3 * 60 * 60)) * 86400 / (2 * pi)) % 86400
                            Ste = Ste * (2 * pi) / 86400

                            Xi = X0 * math.cos(Ste) - Y0 * math.sin(Ste)
                            Yi = X0 * math.sin(Ste) + Y0 * math.cos(Ste)
                            Zi = Z0
                            Xdoti = Xdot * math.cos(Ste) - Ydot * math.sin(Ste) - we * Y0
                            Ydoti = Xdot * math.sin(Ste) + Ydot * math.cos(Ste) + we * X0
                            Zdoti = Zdot

                            w10 = X0
                            w20 = Y0
                            w30 = Z0
                            w40 = Xdot
                            w50 = Ydot
                            w60 = Zdot
                            deltat = i - tb - lps
                            h = sign(deltat)*1 #STEP
                            dts = (TauN - GammaN*(float(i) - float(tb)))
                            aux = 0
                            if abs(deltat)-abs(h)<0:
                                aux = 1
                            for j in range(int(abs(deltat)/1)+ aux):
                                if abs(deltat)-(j+1)*abs(h)< abs(h):
                                    h = sign(h)*(abs(deltat)-j*abs(h))

                                r = math.sqrt((w10 ** 2) + (w20 ** 2) + (w30 ** 2));

                                k11 = h * (w40);
                                k12 = h * (w50);
                                k13 = h * (w60);

                                k14 = h * ((-mu * w10) / (r ** 3) + (
                                        ((3.0 / 2.0) * c20 * mu * (ae ** 2) * w10) / (r ** 5)) * (
                                                1 - 5 * (w30 ** 2) / (r ** 2)) + (
                                                we ** 2) * w10 + 2 * we * w50 + xdotdot)

                                k15 = h * (-mu * w20 / (r ** 3) + (
                                        ((3.0 / 2.0) * c20 * mu * (ae ** 2) * w20) / (r ** 5)) * (
                                                1 - (w30 ** 2) * 5 / (r ** 2)) + (
                                                we ** 2) * w20 - 2 * we * w40 + ydotdot)

                                k16 = h * (-mu * w30 / (r ** 3) + (
                                        ((3 / 2) * c20 * mu * (ae ** 2) * w30) / (r ** 5)) * (
                                                3 - 5 * (w30 ** 2) / (r ** 2)) + zdotdot)

                                ##aqui##
                                k21 = h * (w40 + 0.5 * k14);
                                k22 = h * (w50 + 0.5 * k15);
                                k23 = h * (w60 + 0.5 * k16);

                                k24 = h * ((-mu * (w10 + 0.5 * k11)) / (r ** 3) + (
                                        3 / 2 * c20 * mu * (ae ** 2) * (w10 + 0.5 * k11)) / (r ** 5) * (
                                                1 - 5 * (w30 + 0.5 * k13) ** 2 / (r ** 2)) + (we ** 2) * (
                                                w10 + 0.5 * k11) + 2 * we * (w50 + 0.5 * k15) + xdotdot);
                                # aqui##
                                k25 = h * (-mu / (r ** 3) * (w20 + 0.5 * k12) + 3 / 2 * c20 * mu * (ae ** 2) / (
                                        r ** 5) * (w20 + 0.5 * k12) * (
                                                1 - 5 / (r ** 2) * (w30 + 0.5 * k13) ** 2) + (we ** 2) * (
                                                w20 + 0.5 * k12) - 2 * we * (w40 + 0.5 * k14) + ydotdot);

                                k26 = h * (-mu / (r ** 3) * (w30 + 0.5 * k13) + 3 / 2 * c20 * mu * (ae ** 2) / (
                                        r ** 5) * (w30 + 0.5 * k13) * (
                                                3 - 5 / (r ** 2) * ((w30 + 0.5 * k13) ** 2)) + zdotdot);
                                # aqui#
                                k31 = h * (w40 + 0.5 * k24);
                                k32 = h * (w50 + 0.5 * k25);
                                k33 = h * (w60 + 0.5 * k26);

                                k34 = h * (-mu * (w10 + 0.5 * k21) / (r ** 3) + 3.0 / 2.0 * c20 * mu * (ae ** 2) * (
                                        w10 + 0.5 * k21) / (r ** 5) * (1 - 5 * (w30 + 0.5 * k23) ** 2) / (
                                                r ** 2) + (we ** 2) * (w10 + 0.5 * k21) + 2 * we * (
                                                w50 + 0.5 * k25) + xdotdot)

                                k35 = h * (-mu / (r ** 3) * (w20 + 0.5 * k22) + 3 / 2 * c20 * mu * (ae ** 2) / (
                                        r ** 5) * (w20 + 0.5 * k22) * (
                                                1 - 5 / (r ** 2) * (w30 + 0.5 * k23) ** 2) + (we ** 2) * (
                                                w20 + 0.5 * k22) - 2 * we * (w40 + 0.5 * k24) + ydotdot);

                                k36 = h * (-mu / (r ** 3) * (w30 + 0.5 * k23) + 3 / 2 * c20 * mu * (ae ** 2) / (
                                        r ** 5) * (w30 + 0.5 * k23) * (
                                                3 - 5 / (r ** 2) * ((w30 + 0.5 * k23) ** 2)) + zdotdot);

                                k41 = h * (w40 + k34);
                                k42 = h * (w50 + k35);
                                k43 = h * (w60 + k36);

                                k44 = h * (-mu / (r ** 3) * (w10 + k31) + 3 / 2 * c20 * mu * (ae ** 2) / (r ** 5) * (
                                        w10 + k31) * (1 - 5 / (r ** 2) * ((w30 + k33) ** 2)) + (we ** 2) * (
                                                w10 + k31) + 2 * we * (w50 + k35) + xdotdot)

                                k45 = h * (-mu / (r ** 3) * (w20 + k32) + 3 / 2 * c20 * mu * (ae ** 2) / (r ** 5) * (
                                        w20 + k32) * (1 - 5 / (r ** 2) * ((w30 + k33) ** 2)) + (we ** 2) * (
                                                w20 + k32) - 2 * we * (w40 + k34) + ydotdot)

                                k46 = h * (-mu / (r ** 3) * (w30 + k33) + 3 / 2 * c20 * mu * (ae ** 2) / (r ** 5) * (
                                        w30 + k33) * (3 - 5 / (r ** 2) * ((w30 + k33) ** 2)) + zdotdot);

                                w11 = w10 + (k11 + 2 * k21 + 2 * k31 + k41) / 6;
                                w21 = w20 + (k12 + 2 * k22 + 2 * k32 + k42) / 6;
                                w31 = w30 + (k13 + 2 * k23 + 2 * k33 + k43) / 6;
                                w41 = w40 + (k14 + 2 * k24 + 2 * k34 + k44) / 6;
                                w51 = w50 + (k15 + 2 * k25 + 2 * k35 + k45) / 6;
                                w61 = w60 + (k16 + 2 * k26 + 2 * k36 + k46) / 6;

                                w10 = w11;
                                w20 = w21;
                                w30 = w31;
                                w40 = w41;
                                w50 = w51;
                                w60 = w61;
                            vx = w40 * 1000
                            vy = w50 * 1000
                            vz = w60 * 1000
                            x = w10 * 1000 - 0.36
                            y = w20 * 1000 - 0.08
                            z = w30 * 1000 - 0.18

                            ml = gl
                            seconds = i % 86400+900
                            seconds = seconds % (24 * 3600)
                            hour = seconds // 3600
                            seconds %= 3600
                            minutes = seconds // 60
                            seconds %= 60
                            for l in range(len(ml)):
                                if ml[l].startswith('*'):
                                    c = ml[l].split()
                                    if c[0] == 'EOF':
                                        break
                                    if ((year == int(c[1])) and (month == int(c[2])) and (day == int(c[3])) and (hour == int(c[4])) and (minutes == int(c[5])) and (seconds == int(float((c[6]))))):
                                        m = l;
                                        s1 = ''
                                        while s1 != '*' and s != 'EOF':
                                            cor1 = ml[m + 1].split()
                                            s1 = cor1[0]
                                            if sat in s1:
                                                xp = float(cor1[1]) * 1000
                                                yp = float(cor1[2]) * 1000
                                                zp = float(cor1[3]) * 1000
                                                dtp = float(cor1[4])
                                            m = m + 1

                            V3 = [float(vx), float(vy), float(vz)]
                            V2 = [float(x), float(y), float(z)]
                            V2m = linalg.norm(V2)

                            ru = V2 / V2m  # #Radial component
                            nu = cross(V2, V3) / linalg.norm(cross(V2, V3))  # Along track
                            tu = cross(nu, ru)  # Cross-track

                            V2 = [[xp], [yp], [zp]]  # Precise coordinates
                            V1 = [[x], [y], [z]]  # broadcast coordinate
                            Rct = matrix([ru, tu, nu])  # rotation matrix for RTN
                            Diff2 = Rct * V2 - Rct * V1
                            dradial2 = float(Diff2[0])
                            dalong2 = float(Diff2[1])
                            dcross2 = float(Diff2[2])
                            mod = '{:15.2f}'.format(math.sqrt(dradial2 ** 2 + dalong2 ** 2 + dcross2 ** 2))
                            dradial2 = '{:15.2f}'.format(dradial2)
                            dalong2 = '{:15.2f}'.format(dalong2)
                            dcross2 = '{:15.2f}'.format(dcross2)
                            deltat = '{:10.10f}'.format((dts - dtp * 10 ** -6) * 299792458)
                            cor.write('P' + sat + ' , ' + str(i+900) + ',' + str(dradial2) + ',' + str(dalong2) + ',' + str(dcross2) + ','+ str(mod)+ ',' + str(deltat) +','+ str(int(health))+'\n')