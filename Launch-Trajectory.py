# trajectory model

#  units: SI
#           pos:   altitude           m
#           vel:   velocity           m/s
#           a:     speed of sound     m/s
#           acc:   acceleration       m/s^2
#           g0:    gravity            m/s^2
#           u:     kinem. viscosity   m^2/s
#           T:     temperature        °K
#           P:     pressure           Pa
#           m:     mass               kg
#           s:     frontal area       m^2
#           rho:   density            kg/m^3
#           R:     gas constant       N-m/kg-K
#           d:     drag               N
#           cD:    coef of drag       unitless

# TODO: coefficient of drag function
# TODO: thrust vector
# TODO: drag vector

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# constants
gamma = 1.4
R = 287.05287               # gas constant, N-m/kg-K
g0 = 9.80665                # m/s^2
RE = 6356766                # radius of the Earth, m
Bs = 1.458e-6               # N-s/m2 K1/2
S = 110.4                   # K
azimuth = 45                # degrees
angleOfAttack = 45          # degrees
theta = 90 - angleOfAttack  # degrees
phi = 90 - azimuth          # degrees

# constant performance parameters
coeffDrag = np.random.uniform(.1, .7)
s = np.random.uniform(.1, .7)
burnTime = 1.8  # s

# initial values
initialPosX = 0    # m
initialPosY = 0    # m
initialPosZ = 210  # m
initialVelX = .001*np.sin(theta*(np.pi/180))*np.cos(phi*(np.pi/180))  # m/s
initialVelY = .001*np.sin(theta*(np.pi/180))*np.sin(phi*(np.pi/180))  # m/s
initialVelZ = .001*np.cos(theta*(np.pi/180))                          # m/s


def atmosphericConditions(altitude):
    # calculate gas properties in earth's atmosphere
    # setup indices for each atmospheric layer
    #                          index    lapse rate   base Temp        alt                  base pressure
    #                          i        Ki (°C/m)    Ti (°K)          Hi (m)               P (Pa)
    atmLayerTable = np.array([[1,       -.0065,      288.15,          0,                   101325],
                              [2,       0,           216.65,          11000,               22632.0400950078],
                              [3,      .001,         216.65,          20000,               5474.87742428105],
                              [4,       .0028,       228.65,          32000,               868.015776620216],
                              [5,       0,           270.65,          47000,               110.90577336731],
                              [6,       -.0028,      270.65,          51000,               66.9385281211797],
                              [7,       -.002,       214.65,          71000,               3.9563921603966],
                              [8,       0,           186.94590831019, 84852.0458449057,    0.373377173762337]])

    atmLayerK = atmLayerTable[:, 1]  # °K/m
    atmLayerT = atmLayerTable[:, 2]	 # °K
    atmLayerH = atmLayerTable[:, 3]	 # m
    atmLayerP = atmLayerTable[:, 4]  # Pa

    altitudeMax = 90000  # m

    # initialize P and T
    P = 0
    T = 0

    # troposphere
    if altitude <= atmLayerH[1]:
        i = 0
        TonTi = 1 + atmLayerK[i] * (altitude - atmLayerH[i]) / atmLayerT[i]
        T = TonTi * atmLayerT[i]
        PonPi = TonTi ** (-g0 / (atmLayerK[i] * R))
        P = atmLayerP[i] * PonPi

    # tropopause
    if (altitude <= atmLayerH[2]) & (altitude > atmLayerH[1]):
        i = 1
        T = atmLayerT[i]
        PonPi = np.exp(-g0 * (altitude - atmLayerH[i]) / (atmLayerT[i] * R))
        P = PonPi * atmLayerP[i]

    # stratosphere 1
    if (altitude <= atmLayerH[3]) & (altitude > atmLayerH[2]):
        i = 2
        TonTi = 1 + atmLayerK[i] * (altitude - atmLayerH[i]) / atmLayerT[i]
        T = TonTi*atmLayerT[i]
        PonPi = TonTi ** (-g0 / (atmLayerK[i] * R))
        P = PonPi * atmLayerP[i]

    # stratosphere 2
    if (altitude <= atmLayerH[4]) & (altitude > atmLayerH[3]):
        i = 3
        TonTi = 1 + atmLayerK[i] * (altitude - atmLayerH[i]) / atmLayerT[i]
        T = TonTi*atmLayerT[i]
        PonPi = TonTi ** (-g0 / (atmLayerK[i] * R))
        P = PonPi * atmLayerP[i]

    # stratopause
    if (altitude <= atmLayerH[5]) & (altitude > atmLayerH[4]):
        i = 4
        T = atmLayerT[i]
        PonPi = np.exp(-g0 * (altitude - atmLayerH[i]) / (atmLayerT[i] * R))
        P = PonPi * atmLayerP[i]

    # mesosphere 1
    if (altitude <= atmLayerH[6]) & (altitude > atmLayerH[5]):
        i = 5
        TonTi = 1 + atmLayerK[i] * (altitude - atmLayerH[i]) / atmLayerT[i]
        T = TonTi*atmLayerT[i]
        PonPi = TonTi ** (-g0 / (atmLayerK[i] * R))
        P = PonPi * atmLayerP[i]

    # mesosphere 2
    if (altitude <= atmLayerH[7]) & (altitude > atmLayerH[6]):
        i = 6
        TonTi = 1 + atmLayerK[i] * (altitude - atmLayerH[i]) / atmLayerT[i]
        T = TonTi*atmLayerT[i]
        PonPi = TonTi ** (-g0 / (atmLayerK[i] * R))
        P = PonPi * atmLayerP[i]

    # mesopause
    if (altitude <= altitudeMax) & (altitude > atmLayerH[7]):
        i = 7
        T = atmLayerT[i]
        PonPi = np.exp(-g0 * (altitude - atmLayerH[i]) / (atmLayerT[i] * R))
        P = PonPi * atmLayerP[i]

    # thermosphere
    if altitude > altitudeMax:
        print('WARNING: altitude above atmospheric upper limit')
        T = 0
        P = 0

    # calc density (rho), acoustic speed (a), and kinematic viscosity (u)
    rho = P / (T * R)                     # kg/m^3
    a = (gamma * R * T) ** 0.5            # m/s
    u = (Bs * (T**1.5) / (T + S)) / rho   # m^2/s
    return T, P, rho, a, u


def calcDrag(s, rho, velocity, cD):
    drag = 0.5 * rho * cD * s * velocity**2  # N
    return drag


def accX(thrust, drag, mass):
    aX = (thrust - drag) / mass
    return aX


def accY(thrust, drag, mass):
    aY = (thrust - drag) / mass
    return aY


def accZ(thrust, drag, mass, gravity):
    aZ = (thrust - drag) / mass - gravity
    return aZ


def odeModelX(y, t):
    x, vX = y
    dxdt = [vX, aX]
    return dxdt


def odeModelY(y, t):
    x, vY = y
    dydt = [vY, aY]
    return dydt


def odeModelZ(y, t):
    x, vZ = y
    dzdt = [vZ, aZ]
    return dzdt


def calcMass(time):
    massProp = 820
    massInert = 65.6
    burnRate = np.absolute((massProp - massInert) / burnTime)  # kg/s

    mass = massProp + massInert
    if mass < massInert:
        mass = massProp + massInert - burnRate * time
    else:
        mass = massInert
    return mass  # kg


def plotting(time, posX, posZ):
    plt.plot(time, posX, label='x(t)')
    plt.plot(time, posZ, label='z(t)')
    plt.legend(loc='best')
    plt.xlabel('time (s)')
    plt.ylabel('distance (m)')
    plt.grid()
    plt.annotate(f'{posX[-1]:.2f} m', (time[-1], posX[-1]), textcoords="offset points", xytext=(0, 0), ha='right')
    plt.show()


if __name__ == '__main__':
    t = [0]
    time = [0]
    j = 0
    dt = .01

    x0 = [initialPosX, initialVelX]
    y0 = [initialPosY, initialVelY]
    z0 = [initialPosZ, initialVelZ]

    posX = [initialPosX]
    posY = [initialPosY]
    posZ = [initialPosZ]

    velX = [initialVelX]
    velY = [initialVelY]
    velZ = [initialVelZ]

    while posZ[-1] > 0:
        if j < 1.8:
            thrust = 215250
            t = [j, j + dt]
            temp, press, density, acousticSpeed, kinVisc = atmosphericConditions(z0[0])

            # drag
            dragX = calcDrag(s, density, x0[1], coeffDrag)
            dragY = calcDrag(s, density, y0[1], coeffDrag)
            dragZ = calcDrag(s, density, z0[1], coeffDrag)

            # acceleration
            aX = accX(thrust, dragX, calcMass(t[1]))
            aY = accY(thrust, dragY, calcMass(t[1]))
            aZ = accZ(thrust, dragZ, calcMass(t[1]), g0)

            # integration
            solX = odeint(odeModelX, x0, t)
            solY = odeint(odeModelY, y0, t)
            solZ = odeint(odeModelZ, z0, t)

            # update initial values
            x0 = [solX[1, 0], solX[1, 1]]
            y0 = [solY[1, 0], solY[1, 1]]
            z0 = [solZ[1, 0], solZ[1, 1]]

            # output
            time = np.append([time], t[1])
            posX = np.append([posX], solX[1, 0])
            posY = np.append([posY], solY[1, 0])
            posZ = np.append([posZ], solZ[1, 0])
            velX = np.append([velX], solX[1, 1])
            velY = np.append([velY], solY[1, 1])
            velZ = np.append([velZ], solZ[1, 1])

            j = j + dt

        else:
            thrust = 0
            t = [j, j + dt]
            temp, press, density, acousticSpeed, kinVisc = atmosphericConditions(z0[0])

            # drag
            dragX = calcDrag(s, density, x0[1], coeffDrag)
            dragY = calcDrag(s, density, y0[1], coeffDrag)
            dragZ = calcDrag(s, density, z0[1], coeffDrag)

            # acceleration
            aX = accX(thrust, dragX, calcMass(t[1]))
            aY = accY(thrust, dragY, calcMass(t[1]))
            aZ = accZ(thrust, dragZ, calcMass(t[1]), g0)

            # integration
            solX = odeint(odeModelX, x0, t)
            solY = odeint(odeModelY, y0, t)
            solZ = odeint(odeModelZ, z0, t)

            # update initial values
            x0 = [solX[1, 0], solX[1, 1]]
            y0 = [solY[1, 0], solY[1, 1]]
            z0 = [solZ[1, 0], solZ[1, 1]]

            # output
            time = np.append([time], t[1])
            posX = np.append([posX], solX[1, 0])
            posY = np.append([posY], solY[1, 0])
            posZ = np.append([posZ], solZ[1, 0])
            velX = np.append([velX], solX[1, 1])
            velY = np.append([velY], solY[1, 1])
            velZ = np.append([velZ], solZ[1, 1])

            j = j + dt

    outputPosition = np.column_stack((posX, posY, posZ))
    outputVelocity = np.column_stack((velX, velY, velZ))
    plotting(time, posX, posZ)
