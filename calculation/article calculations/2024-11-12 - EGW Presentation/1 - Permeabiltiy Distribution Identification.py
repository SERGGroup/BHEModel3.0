# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_classes.constant import ARTICLE_CALCULATION_DIR
from scipy.stats import norm, lognorm
import matplotlib.pyplot as plt
import numpy as np
import arviz as az
import pymc as pm
import openpyxl
import os


# %%------------   IMPORT DATA                            -----------------------------------------------------------> #
file_path = os.path.join(

    ARTICLE_CALCULATION_DIR, '2024-11-12 - EGW Presentation',
    '0 - Res', 'Dati Gavorrano_HOCLOOP.xlsx'

)

# Load the workbook
workbook = openpyxl.load_workbook(file_path)

# Select the active worksheet (or specify a sheet name)
sheet = workbook['Foglio1']

k_data = list()

for cell in sheet['I']:

    if type(cell.value) is float:

        k_data.append(cell.value)

k_data = np.array(k_data)
k_log_data = np.log10(k_data)


# %%------------   PERFORM BAYESIAN FITTING               -----------------------------------------------------------> #
use_log_data = True

if use_log_data:
    data = k_log_data
    title = "log(Permeability) Distribution"

else:
    data = k_data
    title = "Permeability Distribution"

base_data_model = pm.Model()
with base_data_model as model:

    # Priors for unknown model parameters
    mu = pm.Uniform('mu', lower=np.min(data), upper=np.max(data))
    sigma = pm.Uniform('sigma', lower=0, upper=(np.max(data) - np.min(data)))

    # Likelihood (sampling distribution) of observations
    likelihood_normal = pm.Normal('obs', mu=mu, sigma=sigma, observed=data)

    # Inference: draw posterior samples using MCMC
    idata = pm.sample(draws=3000)


# %%------------   FIT POSTERIOR                          -----------------------------------------------------------> #
fig, axs = plt.subplots(2, 2, figsize=(10, 7))
az.plot_trace(idata, axes=axs)
posterior_samples = idata.posterior
mu_samples = posterior_samples['mu'].values.flatten()
mu, std = norm.fit(mu_samples)
mu_posterior = {

    "mu": mu,
    "std": std

}
x_values = np.linspace(min(mu_samples), max(mu_samples), 1000)
pdf_fitted = norm.pdf(x_values, mu, std)
axs[0, 0].plot(x_values, pdf_fitted, 'r-', lw=2, label='Fitted norm PDF')


sigma_samples = posterior_samples['sigma'].values.flatten()
shape, loc, scale = lognorm.fit(sigma_samples, floc=0)
sigma_posterior = {

    "shape": shape,
    "loc": loc,
    "scale": scale,
    "mode": np.exp(np.log(scale) - shape**2)

}

x_values = np.linspace(min(sigma_samples), max(sigma_samples), 1000)
pdf_fitted = lognorm.pdf(x_values, shape, loc, scale)
axs[1, 0].plot(x_values, pdf_fitted, 'r-', lw=2, label='Fitted lognormal PDF')

plt.suptitle(title, fontsize=18)
plt.tight_layout(pad=1)
plt.show()


# %%------------   SHOW POSSIBLE GAUSSIANS               -----------------------------------------------------------> #
n_curves = 1000

# Extract posterior samples of mu and sigma
posterior_samples = idata.posterior
mu_samples = posterior_samples['mu'].values.flatten()
sigma_samples = posterior_samples['sigma'].values.flatten()

if use_log_data:
    x_values = np.power(10, np.linspace(np.min(data)*0.88, np.max(data)*1.5, 1000))
    x_pdf = np.log10(x_values)

else:
    x_values = np.linspace(0, np.max(data)*1.4, 1000)
    x_pdf = x_values

# Plot the data histogram and the Gaussian curves
plt.figure(figsize=(10, 6))

# Plot the histogram of the data
# plt.hist(data, bins=30, density=True, alpha=0.5, label='Data')

# Plot 10 random Gaussian curves using posterior samples
for i in np.random.randint(0, len(mu_samples), n_curves):

    mu = mu_samples[i]
    sigma = sigma_samples[i]

    # Compute the Gaussian curve
    y_values = norm.pdf(x_pdf, loc=mu, scale=sigma)

    # Plot the curve
    plt.plot(x_values, y_values, alpha=0.01, color="tab:blue")

pdf_fitted= norm.pdf(x_pdf, mu_posterior["mu"], sigma_posterior["mode"])
plt.plot(x_values, pdf_fitted,  lw=3, alpha=1, color="tab:blue")

# Add labels and title
plt.title(f'{title} with Posterior Gaussian Curves')
plt.xlabel('Permeability [m$^2$]')
plt.ylabel('Probability Density')

plt.show()