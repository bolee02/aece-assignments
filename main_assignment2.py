import numpy as np
import matplotlib.pyplot as plt

# constants
eta_0 = 0.3  # unitless
eta_1 = 0.4  # unitless

EI_H2O_kero = 1250.  # g/kg
EI_H2O_H2 = 8940.  # g/kg

M_air = 29.  # kg/(g mol)
M_H2O = 18.
c_p = 1004.  # J/(kg K)

Q_kero = 43.2  # MJ/kg
Q_H2 = 120.  # MJ/kg


temp = -48.4 + 273.15  # K
press = 274 * 10 ** 2  # Pa


a_water = [-6.0969385 * 10 ** 3, 2.12409642 * 10, -2.71193 * 10 ** -2, 1.673952 * 10 ** -5, 2.433502]
a_ice = [-6.0245282 * 10 ** 3, 2.932707 * 10, 1.0613868 * 10 ** -2, -1.3198825 * 10 ** -5, -4.9382577 * 10 ** -1]


def saturation_press(T):
    exponent_water = 0
    exponent_ice = 0
    for i in range(4):
        exponent_water += a_water[i] * T ** (i - 1)
        exponent_ice += a_ice[i] * T ** (i - 1)
    return np.exp(exponent_water + a_water[4] * np.log(T)), np.exp(exponent_ice + a_ice[4] * np.log(T))


temp_list = np.arange(215, 280, dtype=float)
sat_ice = np.array(list(map(saturation_press, temp_list)))[:,1]
sat_water = np.array(list(map(saturation_press, temp_list)))[:,0]

p_a = saturation_press(temp)[1]
G_0 = press/1000 * c_p * M_air/M_H2O * EI_H2O_kero/((1 - eta_0) * Q_kero * 10 ** 6)
G_1 = press/1000 * c_p * M_air/M_H2O * EI_H2O_kero/((1 - eta_1) * Q_kero * 10 ** 6)
G_2 = press/1000 * c_p * M_air/M_H2O * EI_H2O_H2/((1 - eta_0) * Q_H2 * 10 ** 6)
p_0, p_1, p_2 = [], [], []
temp_list_ac = np.arange(round(temp), 280, dtype=float)
for i in range(len(temp_list_ac)):
    p_0.append(p_a + G_0 * (temp_list_ac[i] - temp))
    p_1.append(p_a + G_1 * (temp_list_ac[i] - temp))
    p_2.append(p_a + G_2 * (temp_list_ac[i] - temp))

plt.plot(temp_list, sat_water, color='darkblue')
plt.plot(temp_list, sat_ice, color='lightblue')
plt.plot(temp_list_ac, p_0, color='orange', linestyle='dashed')
plt.plot(temp_list_ac, p_1, color='red', linestyle='dashed')
plt.plot(temp_list_ac, p_2, color='darkred', linestyle='dashed')
plt.xlabel('Temperature (K)')
plt.ylabel('Water Vapour Pressure (Pa)')
plt.legend(['Saturation w.r.t. water', 'Saturation w.r.t. ice',
            'Mixing line - kero, eta = 0.3', 'Mixing line - kero, eta = 0.4',
            'Mixing line - H2, eta = 0.3'])
plt.title("Schmidt-Appleman Plots")
#plt.ylim(0, 200)
plt.grid()
plt.show()
