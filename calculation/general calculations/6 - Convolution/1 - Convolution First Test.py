# %%-------------------------   IMPORT MODULES                      -------------------------------------------------> #
import matplotlib.pyplot as plt
from scipy.special import kn
from tqdm import tqdm
import numpy as np


class ReservoirProperties:

    grad = 0.05      # [°C/m]
    c_rock = 902.67  # [J/(kg K)]
    rho_rock = 2600  # [kg/m^3]
    k_rock = 2.423   # [W/(m K)]
    pe = 1.  # [-]

    __corr = np.array([

        [0.6085, 0.2760, 0.2165],
        [1.3465, 0.4151, 0.7052],
        [0.3777, 0.2792, 0.2195],
        [0.3324, 2.8015, 2.6310],
        [0.2082, 0.0962, 0.5448]

    ])

    @property
    def alpha_rock(self) -> float:
        #   alpha_rocks -> rock thermal diffusivity in [m^2/s]
        return self.k_rock / (self.rho_rock * self.c_rock)

    @property
    def grad_km(self) -> float:
        return self.grad * 1E3

    @grad_km.setter
    def grad_km(self, grad_km: float):
        self.grad = grad_km / 1E3

    def evaluate_rel_resistance(self, times: [float], d: float) -> [float]:

        times_nd = self.get_nd_times(times, d)
        f = self.evaluate_f(times_nd)

        return d / (2 * self.k_rock * f)

    def evaluate_f(self, times_nd: [float]) -> [float]:

        params = self.__corr[:, 1] * np.power(self.pe, self.__corr[:, 0]) + self.__corr[:, 2]

        f = np.empty(times_nd.shape)
        f[:] = params[0]

        for i in range(int(np.floor((len(params) - 1) / 2))):
            f += kn(i, params[1 + i * 2] * np.power(times_nd, params[2 + i * 2]))

        return f

    def get_nd_times(self, times: [float], d: float) -> [float]:
        return np.array(times) * (4 * self.alpha_rock) / np.power(d, 2)

    def get_times(self, nd_times: [float], d: float) -> [float]:
        return np.array(nd_times) / (4 * self.alpha_rock) * np.power(d, 2)


# %%-------------------------   DEFINE VARIABLES                    -------------------------------------------------> #
cp_h2o = 4184   # (J / kg °C)
rho_h2o = 1000  # (kg / m^3)
r_well = 0.05   # (m)
l_well = 1      # (m)
T_0 = 20        # (°C)
q_in_val = (80, 20)     # (W)

a_well = 2 * np.pi * r_well * l_well
m_h2o = rho_h2o * (l_well * np.pi * r_well ** 2)
c_h2o = cp_h2o * m_h2o

rp = ReservoirProperties()
rp.k_rock = 2.1     # (W/m °C)
rp.c_rock = 1000    # (J/kg °C)
rp.rho_rock = 2650  # (kg/m^3)
rp.pe = 0           # (-)

times = np.linspace(0, 5E6, 20001)
dt = times[1] - times[0]

q_net = np.zeros((len(times), 3))
dT_in = np.zeros((len(times), 3))
T_in = np.zeros((len(times), 3))

times_floor = np.remainder(np.floor(times / 1E6), 2)
i_even = np.where(times_floor == 0)
i_odd = np.where(times_floor == 1)

q_in = np.ones(times.shape)
q_in[i_even] *= q_in_val[0]
q_in[i_odd] *= q_in_val[1]

base_res = a_well / rp.evaluate_rel_resistance(times, r_well * 2)   # (W/°C)


# %%-------------------------   EVALUATE CONVOLUTION                -------------------------------------------------> #
pbar = tqdm(desc="Calculating ", total=len(times))
for i in range(1, len(times)):

    dT_in[i, :] = q_net[i - 1, :] * dt / c_h2o
    T_in[i, :] = T_in[i - 1, :] + dT_in[i, :]

    # NO HEAT EXCHANGE
    q_net[i, 0] = q_in[i]

    # NO CONVOLUTION
    q_net[i, 1] = q_in[i] - base_res[i] * T_in[i, 1]

    # CONVOLUTION
    q_net[i, 2] = q_in[i]
    for j in range(i - 1):
        k = j + 1
        q_net[i, 2] -= base_res[k] * dT_in[i - k, 2]

    pbar.update(1)

pbar.close()


# %%-------------------------   PLOT TEMPERATURES                   -------------------------------------------------> #
plt.plot(times, T_in[:, 0] + 20, label="No HE")
plt.plot(times, T_in[:, 1] + 20, label="No Convolution")
plt.plot(times, T_in[:, 2] + 20, label="Convolution")

plt.ylim((20, 50))
plt.xlim((100, max(times)))
plt.xscale("log")
plt.legend()
plt.show()


# %%-------------------------   PLOT NET POWER                      -------------------------------------------------> #
plt.plot(times, q_net[:, 0] - 80, label="No HE")
plt.plot(times, q_net[:, 1] - 80, label="No Convolution")
plt.plot(times, q_net[:, 2] - 80, label="Convolution")

plt.xlim((100, max(times)))
plt.xscale("log")
plt.legend()
plt.show()
