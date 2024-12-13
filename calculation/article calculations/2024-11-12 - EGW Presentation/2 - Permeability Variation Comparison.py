# %%-------------------------   IMPORT MODULES                      -------------------------------------------------> #
from scipy.stats import norm, lognorm, skewnorm
from matplotlib import pyplot as plt
from scipy.special import kn
from tqdm import tqdm
import numpy as np
import arviz as az
import pymc as pm


mu_post = {

    'mu': -13.283337196747874,
    'std': 0.08223099935497215

}
sigma_post = {

    'loc': 0,
    'mode': 0.3931599437544563,
    'scale': 0.40219974007567494,
    'shape': 0.15077241274089442

}


class Element:

    t_rock_surf = 10    # [°C]
    grad = 0.05         # [°C/m]
    c_rock = 902.67     # [J/(kg K)]
    rho_rock = 2600     # [kg/m^3]
    k_rock = 2.423      # [W/(m K)]
    pe = 1.             # [-]

    l_segment = 1       # [m]
    depth = 2000        # [m]
    d = 0.15            # [m]
    m_dot = 1           # [kg/s]

    __corr = np.array([

        [0.6085, 0.2760, 0.2165],
        [1.3465, 0.4151, 0.7052],
        [0.3777, 0.2792, 0.2195],
        [0.3324, 2.8015, 2.6310],
        [0.2082, 0.0962, 0.5448]

    ])

    def __init__(self, previous=None):

        self.previous = previous
        self.next = None

        if previous is not None:
            self.previous.next = self

        self.__T_prev = 0.
        self.__T_in = 0.
        self.t_curr = 0.

        self.T_in_list = list()
        self.t_list = list()

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

    @property
    def t_rock_res(self) -> float:

        return self.t_rock_surf + self.grad * self.depth

    @property
    def surf_area(self) -> float:

        return self.l_segment * np.pi * self.d

    def evaluate_rel_resistance(self, times: [float]) -> [float]:

        times_nd = self.get_nd_times(times)
        f = self.evaluate_f(times_nd)

        return self.d / (2 * self.k_rock * f)

    def evaluate_f(self, times_nd: [float]) -> [float]:

        params = self.__corr[:, 1] * np.power(self.pe, self.__corr[:, 0]) + self.__corr[:, 2]

        f = np.empty(times_nd.shape)
        f[:] = params[0]

        for i in range(int(np.floor((len(params) - 1) / 2))):
            f += kn(i, params[1 + i * 2] * np.power(times_nd, params[2 + i * 2]))

        return f

    def get_nd_times(self, times: [float]) -> [float]:
        #   alpha_rocks -> rock thermal diffusivity in [m^2/s]
        return np.array(times) * (4 * self.alpha_rock) / np.power(self.d, 2)

    def get_times(self, nd_times: [float]) -> [float]:
        return np.array(nd_times) / (4 * self.alpha_rock) * np.power(self.d, 2)

    def update_new(self, td: float, upward: bool = False):

        if upward:

            new_point = self.previous

        else:

            new_point = self.next

        if new_point is not None:

            f = self.evaluate_f(np.array([td]))[0]
            r_tot = self.d / 2 / (self.k_rock * f)
            cp_h2o = 4184  # (J / kg °C)

            dT_rock = self.t_rock_res - self.T_in
            q_tot = dT_rock / r_tot * self.surf_area     # [W]
            dh_fluid = q_tot / self.m_dot                # [J/kg]

            new_point.T_in = self.T_in + dh_fluid / cp_h2o

    @property
    def T_in(self):

        return self.__T_in

    @T_in.setter
    def T_in(self, T_in: float):

        self.__T_prev = self.__T_in
        self.__T_in = T_in


# %%-------------------------   INIT ELEMENTS                       -------------------------------------------------> #
n_elements = 1000
mean_pe = 0.5

first_element = Element()
first_element.pe = mean_pe
new_element = first_element

for i in range(n_elements - 1):

    new_element = Element(new_element)
    new_element.pe = mean_pe


# %%-------------------------   EVALUATE STANDARD                   -------------------------------------------------> #
times_list = np.logspace(-2, 2, num=5)

for time in times_list:

    T_list = []
    print(time)
    first_element.T_in = 20  # [°C]
    new_element = first_element

    while new_element is not None:

        T_list.append(new_element.T_in)

        new_element.update_new(time)
        new_element = new_element.next

    plt.plot(T_list, label=str(time))

plt.legend()
plt.show()


# %%-------------------------   EVALUATE MODIFIED                   -------------------------------------------------> #
n_calc = 500
t_in = 20       # [°C]
time = 10000

l_sections = 25
n_sections = int(np.ceil(n_elements / l_sections))

mu_list = norm.rvs(loc=mu_post["mu"], scale=mu_post["std"], size=n_calc)
sigma_list = lognorm.rvs(s=sigma_post["shape"], scale=sigma_post["scale"], size=n_calc)
t_out_list = np.empty(n_calc)
t_out_list[:] = np.nan

t_profile = np.empty(n_elements)
pe_profile = np.empty(n_elements)
k_profile = np.empty(n_elements)

limit_t_profiles = np.empty((n_elements, 2))
limit_pe_profiles = np.empty((n_elements, 2))
limit_k_profiles = np.empty((n_elements, 2))

pbar = tqdm(desc="Calculation", total=n_calc)
for i in range(n_calc):

    perm_list = np.power(10, norm.rvs(loc=mu_list[i], scale=sigma_list[i], size=n_sections))
    rel_perm_list = perm_list / np.mean(perm_list)

    first_element.T_in = t_in
    new_element = first_element

    j = 0
    while new_element is not None:

        k = j // l_sections
        new_element.pe = mean_pe * rel_perm_list[k]

        t_profile[j] = new_element.T_in
        pe_profile[j] = new_element.pe
        k_profile[j] = rel_perm_list[k]

        j += 1

        new_element.update_new(time)

        if new_element.next is None:
            break
        else:
            new_element = new_element.next

    t_out_list[i] = new_element.T_in

    if t_out_list[i] == np.nanmin(t_out_list):

        limit_t_profiles[:, 0] = t_profile
        limit_pe_profiles[:, 0] = pe_profile
        limit_k_profiles[:, 0] = k_profile

    elif t_out_list[i] == np.nanmax(t_out_list):

        limit_t_profiles[:, 1] = t_profile
        limit_pe_profiles[:, 1] = pe_profile
        limit_k_profiles[:, 1] = k_profile

    pbar.update(1)

pbar.close()


# %%------------   PERFORM BAYESIAN FITTING               -----------------------------------------------------------> #
t_rel = (t_out_list - t_in)/(first_element.t_rock_res - t_in)
base_data_model = pm.Model()
with base_data_model as model:

    # Prior distributions for the parameters
    mu = pm.Uniform('mu', lower=np.min(t_rel), upper=np.max(t_rel))
    sigma = pm.Uniform('sigma', lower=0, upper=np.max(t_rel)-np.min(t_rel))
    alpha = pm.Normal('alpha', mu=0, sigma=10)

    # Skew-normal likelihood
    likelihood = pm.SkewNormal('likelihood', mu=mu, sigma=sigma, alpha=alpha, observed=t_rel)

    # Sample from the posterior using MCMC
    idata = pm.sample(draws=3000)


# %%------------   FIT POSTERIOR                          -----------------------------------------------------------> #
fig, axs = plt.subplots(3, 2, figsize=(10, 7))
az.plot_trace(idata, axes=axs)
posterior_samples = idata.posterior
params = ["alpha", "mu", "sigma"]
posterior_params = list()

for i, element in enumerate(params):

    samples = posterior_samples[element].values.flatten()
    mu, std = norm.fit(samples)
    posterior_params.append({

        "mu": mu,
        "std": std

    })

    x_values = np.linspace(min(samples), max(samples), 1000)
    pdf_fitted = norm.pdf(x_values, mu, std)
    axs[i, 0].plot(x_values, pdf_fitted, 'r-', lw=2, label='Fitted norm PDF')

plt.tight_layout()
plt.show()


# %%------------   SHOW POSSIBLE DISTRIBUTION             -----------------------------------------------------------> #
n_curves = 1000

# Extract posterior samples of mu and sigma
posterior_samples = idata.posterior
mu_samples = posterior_samples['mu'].values.flatten()
sigma_samples = posterior_samples['sigma'].values.flatten()
alpha_samples = posterior_samples['alpha'].values.flatten()

x_values = np.linspace(np.min(t_rel), np.max(t_rel), 1000)

plt.hist(t_rel, density=True, bins=100, alpha=0.5)

for i in np.random.randint(0, len(alpha_samples), n_curves):

    # Compute the Gaussian curve
    scale = sigma_samples[i]
    shape = alpha_samples[i]
    loc = mu_samples[i]

    y_values = skewnorm.pdf(x_values, shape, loc, scale)

    # Plot the curve
    plt.plot(x_values, y_values, alpha=0.01, color="tab:orange")

scale = posterior_params[2]["mu"]
shape = posterior_params[0]["mu"]
loc = posterior_params[1]["mu"]
pdf_fitted = skewnorm.pdf(x_values, shape, loc, scale)
plt.plot(x_values, pdf_fitted,  lw=3, alpha=1, color="tab:orange")
plt.xlabel("$T_{rel}$ (-)", fontsize=10)
plt.ylabel("Probability Density (-)", fontsize=10)
plt.show()


# %%------------   PLOT PROFILES                         -----------------------------------------------------------> #
plot_min = False
if plot_min:
    k = 0
    title = "Worst Temperature Increase"

else:
    k = 1
    title = "Best Temperature Increase"

x = np.linspace(0, n_elements-1, n_elements) * new_element.l_segment
y = limit_k_profiles[:, k]
alpha_values = (y - limit_k_profiles.min()) / (limit_k_profiles.max() - limit_k_profiles.min())

plt.plot(x, limit_t_profiles[:, k], lw=3)

# Color the background by varying the alpha values
for i in range(len(x)-1):
    plt.axvspan(x[i], x[i+1], color='tab:orange', alpha=alpha_values[i])

plt.title(title, fontsize=12)
plt.xlabel("Length [m]", fontsize=10)
plt.ylabel("Temperature [°C]", fontsize=10)
plt.xlim((min(x), max(x)))
plt.show()